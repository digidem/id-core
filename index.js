var _ = require('lodash');
var rbush = require('rbush');
var d3 = require('./lib/d3');

var iD = {
    data: {},
    areaKeys: {},
    actions: {}
};

iD.data.deprecated = require('iD/data/deprecated');

var locale = { _current: 'en' };

locale.current = function(_) {
    if (!arguments.length) return locale._current;
    if (locale[_] !== undefined) locale._current = _;
    else if (locale[_.split('-')[0]]) locale._current = _.split('-')[0];
    return locale;
};

function t(s, o, loc) {
    loc = loc || locale._current;

    var path = s.split(".").reverse(),
        rep = locale[loc];

    while (rep !== undefined && path.length) rep = rep[path.pop()];

    if (rep !== undefined) {
        if (o) for (var k in o) rep = rep.replace('{' + k + '}', o[k]);
        return rep;
    }

    if (loc !== 'en') {
        return t(s, o, 'en');
    }

    if (o && 'default' in o) {
        return o['default'];
    }

    var missing = 'Missing ' + loc + ' translation: ' + s;
    if (typeof console !== "undefined") console.error(missing);

    return missing;
}
iD.geo = {};

iD.geo.roundCoords = function(c) {
    return [Math.floor(c[0]), Math.floor(c[1])];
};

iD.geo.interp = function(p1, p2, t) {
    return [p1[0] + (p2[0] - p1[0]) * t,
            p1[1] + (p2[1] - p1[1]) * t];
};

// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
iD.geo.cross = function(o, a, b) {
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
};

// http://jsperf.com/id-dist-optimization
iD.geo.euclideanDistance = function(a, b) {
    var x = a[0] - b[0], y = a[1] - b[1];
    return Math.sqrt((x * x) + (y * y));
};

// using WGS84 polar radius (6356752.314245179 m)
// const = 2 * PI * r / 360
iD.geo.latToMeters = function(dLat) {
    return dLat * 110946.257617;
};

// using WGS84 equatorial radius (6378137.0 m)
// const = 2 * PI * r / 360
iD.geo.lonToMeters = function(dLon, atLat) {
    return Math.abs(atLat) >= 90 ? 0 :
        dLon * 111319.490793 * Math.abs(Math.cos(atLat * (Math.PI/180)));
};

// using WGS84 polar radius (6356752.314245179 m)
// const = 2 * PI * r / 360
iD.geo.metersToLat = function(m) {
    return m / 110946.257617;
};

// using WGS84 equatorial radius (6378137.0 m)
// const = 2 * PI * r / 360
iD.geo.metersToLon = function(m, atLat) {
    return Math.abs(atLat) >= 90 ? 0 :
        m / 111319.490793 / Math.abs(Math.cos(atLat * (Math.PI/180)));
};

// Equirectangular approximation of spherical distances on Earth
iD.geo.sphericalDistance = function(a, b) {
    var x = iD.geo.lonToMeters(a[0] - b[0], (a[1] + b[1]) / 2),
        y = iD.geo.latToMeters(a[1] - b[1]);
    return Math.sqrt((x * x) + (y * y));
};

iD.geo.edgeEqual = function(a, b) {
    return (a[0] === b[0] && a[1] === b[1]) ||
        (a[0] === b[1] && a[1] === b[0]);
};

// Return the counterclockwise angle in the range (-pi, pi)
// between the positive X axis and the line intersecting a and b.
iD.geo.angle = function(a, b, projection) {
    a = projection(a.loc);
    b = projection(b.loc);
    return Math.atan2(b[1] - a[1], b[0] - a[0]);
};

// Choose the edge with the minimal distance from `point` to its orthogonal
// projection onto that edge, if such a projection exists, or the distance to
// the closest vertex on that edge. Returns an object with the `index` of the
// chosen edge, the chosen `loc` on that edge, and the `distance` to to it.
iD.geo.chooseEdge = function(nodes, point, projection) {
    var dist = iD.geo.euclideanDistance,
        points = nodes.map(function(n) { return projection(n.loc); }),
        min = Infinity,
        idx, loc;

    function dot(p, q) {
        return p[0] * q[0] + p[1] * q[1];
    }

    for (var i = 0; i < points.length - 1; i++) {
        var o = points[i],
            s = [points[i + 1][0] - o[0],
                 points[i + 1][1] - o[1]],
            v = [point[0] - o[0],
                 point[1] - o[1]],
            proj = dot(v, s) / dot(s, s),
            p;

        if (proj < 0) {
            p = o;
        } else if (proj > 1) {
            p = points[i + 1];
        } else {
            p = [o[0] + proj * s[0], o[1] + proj * s[1]];
        }

        var d = dist(p, point);
        if (d < min) {
            min = d;
            idx = i + 1;
            loc = projection.invert(p);
        }
    }

    return {
        index: idx,
        distance: min,
        loc: loc
    };
};

// Return the intersection point of 2 line segments.
// From https://github.com/pgkelley4/line-segments-intersect
// This uses the vector cross product approach described below:
//  http://stackoverflow.com/a/565282/786339
iD.geo.lineIntersection = function(a, b) {
    function subtractPoints(point1, point2) {
        return [point1[0] - point2[0], point1[1] - point2[1]];
    }
    function crossProduct(point1, point2) {
        return point1[0] * point2[1] - point1[1] * point2[0];
    }

    var p = [a[0][0], a[0][1]],
        p2 = [a[1][0], a[1][1]],
        q = [b[0][0], b[0][1]],
        q2 = [b[1][0], b[1][1]],
        r = subtractPoints(p2, p),
        s = subtractPoints(q2, q),
        uNumerator = crossProduct(subtractPoints(q, p), r),
        denominator = crossProduct(r, s);

    if (uNumerator && denominator) {
        var u = uNumerator / denominator,
            t = crossProduct(subtractPoints(q, p), s) / denominator;

        if ((t >= 0) && (t <= 1) && (u >= 0) && (u <= 1)) {
            return iD.geo.interp(p, p2, t);
        }
    }

    return null;
};

iD.geo.pathIntersections = function(path1, path2) {
    var intersections = [];
    for (var i = 0; i < path1.length - 1; i++) {
        for (var j = 0; j < path2.length - 1; j++) {
            var a = [ path1[i], path1[i+1] ],
                b = [ path2[j], path2[j+1] ],
                hit = iD.geo.lineIntersection(a, b);
            if (hit) intersections.push(hit);
        }
    }
    return intersections;
};

// Return whether point is contained in polygon.
//
// `point` should be a 2-item array of coordinates.
// `polygon` should be an array of 2-item arrays of coordinates.
//
// From https://github.com/substack/point-in-polygon.
// ray-casting algorithm based on
// http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
//
iD.geo.pointInPolygon = function(point, polygon) {
    var x = point[0],
        y = point[1],
        inside = false;

    for (var i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
        var xi = polygon[i][0], yi = polygon[i][1];
        var xj = polygon[j][0], yj = polygon[j][1];

        var intersect = ((yi > y) !== (yj > y)) &&
            (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
        if (intersect) inside = !inside;
    }

    return inside;
};

iD.geo.polygonContainsPolygon = function(outer, inner) {
    return _.every(inner, function(point) {
        return iD.geo.pointInPolygon(point, outer);
    });
};

iD.geo.polygonIntersectsPolygon = function(outer, inner, checkSegments) {
    function testSegments(outer, inner) {
        for (var i = 0; i < outer.length - 1; i++) {
            for (var j = 0; j < inner.length - 1; j++) {
                var a = [ outer[i], outer[i+1] ],
                    b = [ inner[j], inner[j+1] ];
                if (iD.geo.lineIntersection(a, b)) return true;
            }
        }
        return false;
    }

    function testPoints(outer, inner) {
        return _.some(inner, function(point) {
            return iD.geo.pointInPolygon(point, outer);
        });
    }

   return testPoints(outer, inner) || (!!checkSegments && testSegments(outer, inner));
};

iD.geo.pathLength = function(path) {
    var length = 0,
        dx, dy;
    for (var i = 0; i < path.length - 1; i++) {
        dx = path[i][0] - path[i + 1][0];
        dy = path[i][1] - path[i + 1][1];
        length += Math.sqrt(dx * dx + dy * dy);
    }
    return length;
};
iD.geo.Extent = function geoExtent(min, max) {
    if (!(this instanceof iD.geo.Extent)) return new iD.geo.Extent(min, max);
    if (min instanceof iD.geo.Extent) {
        return min;
    } else if (min && min.length === 2 && min[0].length === 2 && min[1].length === 2) {
        this[0] = min[0];
        this[1] = min[1];
    } else {
        this[0] = min        || [ Infinity,  Infinity];
        this[1] = max || min || [-Infinity, -Infinity];
    }
};

iD.geo.Extent.prototype = new Array(2);

_.extend(iD.geo.Extent.prototype, {
    equals: function (obj) {
        return this[0][0] === obj[0][0] &&
            this[0][1] === obj[0][1] &&
            this[1][0] === obj[1][0] &&
            this[1][1] === obj[1][1];
    },

    extend: function(obj) {
        if (!(obj instanceof iD.geo.Extent)) obj = new iD.geo.Extent(obj);
        return iD.geo.Extent([Math.min(obj[0][0], this[0][0]),
                              Math.min(obj[0][1], this[0][1])],
                             [Math.max(obj[1][0], this[1][0]),
                              Math.max(obj[1][1], this[1][1])]);
    },

    _extend: function(extent) {
        this[0][0] = Math.min(extent[0][0], this[0][0]);
        this[0][1] = Math.min(extent[0][1], this[0][1]);
        this[1][0] = Math.max(extent[1][0], this[1][0]);
        this[1][1] = Math.max(extent[1][1], this[1][1]);
    },

    area: function() {
        return Math.abs((this[1][0] - this[0][0]) * (this[1][1] - this[0][1]));
    },

    center: function() {
        return [(this[0][0] + this[1][0]) / 2,
                (this[0][1] + this[1][1]) / 2];
    },

    polygon: function() {
        return [
            [this[0][0], this[0][1]],
            [this[0][0], this[1][1]],
            [this[1][0], this[1][1]],
            [this[1][0], this[0][1]],
            [this[0][0], this[0][1]]
        ];
    },

    contains: function(obj) {
        if (!(obj instanceof iD.geo.Extent)) obj = new iD.geo.Extent(obj);
        return obj[0][0] >= this[0][0] &&
               obj[0][1] >= this[0][1] &&
               obj[1][0] <= this[1][0] &&
               obj[1][1] <= this[1][1];
    },

    intersects: function(obj) {
        if (!(obj instanceof iD.geo.Extent)) obj = new iD.geo.Extent(obj);
        return obj[0][0] <= this[1][0] &&
               obj[0][1] <= this[1][1] &&
               obj[1][0] >= this[0][0] &&
               obj[1][1] >= this[0][1];
    },

    intersection: function(obj) {
        if (!this.intersects(obj)) return new iD.geo.Extent();
        return new iD.geo.Extent([Math.max(obj[0][0], this[0][0]),
                                  Math.max(obj[0][1], this[0][1])],
                                 [Math.min(obj[1][0], this[1][0]),
                                  Math.min(obj[1][1], this[1][1])]);
    },

    percentContainedIn: function(obj) {
        if (!(obj instanceof iD.geo.Extent)) obj = new iD.geo.Extent(obj);
        var a1 = this.intersection(obj).area(),
            a2 = this.area();

        if (a1 === Infinity || a2 === Infinity || a1 === 0 || a2 === 0) {
            return 0;
        } else {
            return a1 / a2;
        }
    },

    padByMeters: function(meters) {
        var dLat = iD.geo.metersToLat(meters),
            dLon = iD.geo.metersToLon(meters, this.center()[1]);
        return iD.geo.Extent(
                [this[0][0] - dLon, this[0][1] - dLat],
                [this[1][0] + dLon, this[1][1] + dLat]);
    },

    toParam: function() {
        return [this[0][0], this[0][1], this[1][0], this[1][1]].join(',');
    }

});
iD.geo.Turn = function(turn) {
    if (!(this instanceof iD.geo.Turn))
        return new iD.geo.Turn(turn);
    _.extend(this, turn);
};

iD.geo.Intersection = function(graph, vertexId) {
    var vertex = graph.entity(vertexId),
        highways = [];

    // Pre-split ways that would need to be split in
    // order to add a restriction. The real split will
    // happen when the restriction is added.
    graph.parentWays(vertex).forEach(function(way) {
        if (!way.tags.highway || way.isArea() || way.isDegenerate())
            return;

        if (way.affix(vertexId)) {
            highways.push(way);
        } else {
            var idx = _.indexOf(way.nodes, vertex.id, 1),
                wayA = iD.Way({id: way.id + '-a', tags: way.tags, nodes: way.nodes.slice(0, idx + 1)}),
                wayB = iD.Way({id: way.id + '-b', tags: way.tags, nodes: way.nodes.slice(idx)});

            graph = graph.replace(wayA);
            graph = graph.replace(wayB);

            highways.push(wayA);
            highways.push(wayB);
        }
    });

    var intersection = {
        highways: highways,
        graph: graph
    };

    intersection.turns = function(fromNodeID) {
        if (!fromNodeID)
            return [];

        var way = _.find(highways, function(way) { return way.contains(fromNodeID); });
        if (way.first() === vertex.id && way.tags.oneway === 'yes')
            return [];
        if (way.last() === vertex.id && way.tags.oneway === '-1')
            return [];

        function withRestriction(turn) {
            graph.parentRelations(graph.entity(turn.from.way)).forEach(function(relation) {
                if (relation.tags.type !== 'restriction')
                    return;

                var f = relation.memberByRole('from'),
                    t = relation.memberByRole('to'),
                    v = relation.memberByRole('via');

                if (f && f.id === turn.from.way &&
                    v && v.id === turn.via.node &&
                    t && t.id === turn.to.way) {
                    turn.restriction = relation.id;
                } else if (/^only_/.test(relation.tags.restriction) &&
                    f && f.id === turn.from.way &&
                    v && v.id === turn.via.node &&
                    t && t.id !== turn.to.way) {
                    turn.restriction = relation.id;
                    turn.indirect_restriction = true;
                }
            });

            return iD.geo.Turn(turn);
        }

        var from = {
                node: way.nodes[way.first() === vertex.id ? 1 : way.nodes.length - 2],
                way: way.id.split(/-(a|b)/)[0]
            },
            via = {node: vertex.id},
            turns = [];

        highways.forEach(function(parent) {
            if (parent === way)
                return;

            var index = parent.nodes.indexOf(vertex.id);

            // backward
            if (parent.first() !== vertex.id && parent.tags.oneway !== 'yes') {
                turns.push(withRestriction({
                    from: from,
                    via: via,
                    to: {node: parent.nodes[index - 1], way: parent.id.split(/-(a|b)/)[0]}
                }));
            }

            // forward
            if (parent.last() !== vertex.id && parent.tags.oneway !== '-1') {
                turns.push(withRestriction({
                    from: from,
                    via: via,
                    to: {node: parent.nodes[index + 1], way: parent.id.split(/-(a|b)/)[0]}
                }));
            }
        });

        // U-turn
        if (way.tags.oneway !== 'yes' && way.tags.oneway !== '-1') {
            turns.push(withRestriction({
                from: from,
                via: via,
                to: from,
                u: true
            }));
        }

        return turns;
    };

    return intersection;
};


iD.geo.inferRestriction = function(graph, from, via, to, projection) {
    var fromWay = graph.entity(from.way),
        fromNode = graph.entity(from.node),
        toWay = graph.entity(to.way),
        toNode = graph.entity(to.node),
        viaNode = graph.entity(via.node),
        fromOneWay = (fromWay.tags.oneway === 'yes' && fromWay.last() === via.node) ||
            (fromWay.tags.oneway === '-1' && fromWay.first() === via.node),
        toOneWay = (toWay.tags.oneway === 'yes' && toWay.first() === via.node) ||
            (toWay.tags.oneway === '-1' && toWay.last() === via.node),
        angle = iD.geo.angle(viaNode, fromNode, projection) -
                iD.geo.angle(viaNode, toNode, projection);

    angle = angle * 180 / Math.PI;

    while (angle < 0)
        angle += 360;

    if (fromNode === toNode)
        return 'no_u_turn';
    if ((angle < 23 || angle > 336) && fromOneWay && toOneWay)
        return 'no_u_turn';
    if (angle < 158)
        return 'no_right_turn';
    if (angle > 202)
        return 'no_left_turn';

    return 'no_straight_on';
};
// For fixing up rendering of multipolygons with tags on the outer member.
// https://github.com/openstreetmap/iD/issues/613
iD.geo.isSimpleMultipolygonOuterMember = function(entity, graph) {
    if (entity.type !== 'way')
        return false;

    var parents = graph.parentRelations(entity);
    if (parents.length !== 1)
        return false;

    var parent = parents[0];
    if (!parent.isMultipolygon() || Object.keys(parent.tags).length > 1)
        return false;

    var members = parent.members, member;
    for (var i = 0; i < members.length; i++) {
        member = members[i];
        if (member.id === entity.id && member.role && member.role !== 'outer')
            return false; // Not outer member
        if (member.id !== entity.id && (!member.role || member.role === 'outer'))
            return false; // Not a simple multipolygon
    }

    return parent;
};

iD.geo.simpleMultipolygonOuterMember = function(entity, graph) {
    if (entity.type !== 'way')
        return false;

    var parents = graph.parentRelations(entity);
    if (parents.length !== 1)
        return false;

    var parent = parents[0];
    if (!parent.isMultipolygon() || Object.keys(parent.tags).length > 1)
        return false;

    var members = parent.members, member, outerMember;
    for (var i = 0; i < members.length; i++) {
        member = members[i];
        if (!member.role || member.role === 'outer') {
            if (outerMember)
                return false; // Not a simple multipolygon
            outerMember = member;
        }
    }

    return outerMember && graph.hasEntity(outerMember.id);
};

// Join `array` into sequences of connecting ways.
//
// Segments which share identical start/end nodes will, as much as possible,
// be connected with each other.
//
// The return value is a nested array. Each constituent array contains elements
// of `array` which have been determined to connect. Each consitituent array
// also has a `nodes` property whose value is an ordered array of member nodes,
// with appropriate order reversal and start/end coordinate de-duplication.
//
// Members of `array` must have, at minimum, `type` and `id` properties.
// Thus either an array of `iD.Way`s or a relation member array may be
// used.
//
// If an member has a `tags` property, its tags will be reversed via
// `iD.actions.Reverse` in the output.
//
// Incomplete members (those for which `graph.hasEntity(element.id)` returns
// false) and non-way members are ignored.
//
iD.geo.joinWays = function(array, graph) {
    var joined = [], member, current, nodes, first, last, i, how, what;

    array = array.filter(function(member) {
        return member.type === 'way' && graph.hasEntity(member.id);
    });

    function resolve(member) {
        return graph.childNodes(graph.entity(member.id));
    }

    function reverse(member) {
        return member.tags ? iD.actions.Reverse(member.id)(graph).entity(member.id) : member;
    }

    while (array.length) {
        member = array.shift();
        current = [member];
        current.nodes = nodes = resolve(member).slice();
        joined.push(current);

        while (array.length && _.first(nodes) !== _.last(nodes)) {
            first = _.first(nodes);
            last  = _.last(nodes);

            for (i = 0; i < array.length; i++) {
                member = array[i];
                what = resolve(member);

                if (last === _.first(what)) {
                    how  = nodes.push;
                    what = what.slice(1);
                    break;
                } else if (last === _.last(what)) {
                    how  = nodes.push;
                    what = what.slice(0, -1).reverse();
                    member = reverse(member);
                    break;
                } else if (first === _.last(what)) {
                    how  = nodes.unshift;
                    what = what.slice(0, -1);
                    break;
                } else if (first === _.first(what)) {
                    how  = nodes.unshift;
                    what = what.slice(1).reverse();
                    member = reverse(member);
                    break;
                } else {
                    what = how = null;
                }
            }

            if (!what)
                break; // No more joinable ways.

            how.apply(current, [member]);
            how.apply(nodes, what);

            array.splice(i, 1);
        }
    }

    return joined;
};
/*
    Bypasses features of D3's default projection stream pipeline that are unnecessary:
    * Antimeridian clipping
    * Spherical rotation
    * Resampling
*/
iD.geo.RawMercator = function () {
    var project = d3.geo.mercator.raw,
        k = 512 / Math.PI, // scale
        x = 0, y = 0, // translate
        clipExtent = [[0, 0], [0, 0]];

    function projection(point) {
        point = project(point[0] * Math.PI / 180, point[1] * Math.PI / 180);
        return [point[0] * k + x, y - point[1] * k];
    }

    projection.invert = function(point) {
        point = project.invert((point[0] - x) / k, (y - point[1]) / k);
        return point && [point[0] * 180 / Math.PI, point[1] * 180 / Math.PI];
    };

    projection.scale = function(_) {
        if (!arguments.length) return k;
        k = +_;
        return projection;
    };

    projection.translate = function(_) {
        if (!arguments.length) return [x, y];
        x = +_[0];
        y = +_[1];
        return projection;
    };

    projection.clipExtent = function(_) {
        if (!arguments.length) return clipExtent;
        clipExtent = _;
        return projection;
    };

    projection.stream = d3.geo.transform({
        point: function(x, y) {
            x = projection([x, y]);
            this.stream.point(x[0], x[1]);
        }
    }).stream;

    return projection;
};
/*
  Order the nodes of a way in reverse order and reverse any direction dependent tags
  other than `oneway`. (We assume that correcting a backwards oneway is the primary
  reason for reversing a way.)

  The following transforms are performed:

    Keys:
          *:right=* ⟺ *:left=*
        *:forward=* ⟺ *:backward=*
       direction=up ⟺ direction=down
         incline=up ⟺ incline=down
            *=right ⟺ *=left

    Relation members:
       role=forward ⟺ role=backward
         role=north ⟺ role=south
          role=east ⟺ role=west

   In addition, numeric-valued `incline` tags are negated.

   The JOSM implementation was used as a guide, but transformations that were of unclear benefit
   or adjusted tags that don't seem to be used in practice were omitted.

   References:
      http://wiki.openstreetmap.org/wiki/Forward_%26_backward,_left_%26_right
      http://wiki.openstreetmap.org/wiki/Key:direction#Steps
      http://wiki.openstreetmap.org/wiki/Key:incline
      http://wiki.openstreetmap.org/wiki/Route#Members
      http://josm.openstreetmap.de/browser/josm/trunk/src/org/openstreetmap/josm/corrector/ReverseWayTagCorrector.java
 */
iD.actions.Reverse = function(wayId) {
    var replacements = [
            [/:right$/, ':left'], [/:left$/, ':right'],
            [/:forward$/, ':backward'], [/:backward$/, ':forward']
        ],
        numeric = /^([+\-]?)(?=[\d.])/,
        roleReversals = {
            forward: 'backward',
            backward: 'forward',
            north: 'south',
            south: 'north',
            east: 'west',
            west: 'east'
        };

    function reverseKey(key) {
        for (var i = 0; i < replacements.length; ++i) {
            var replacement = replacements[i];
            if (replacement[0].test(key)) {
                return key.replace(replacement[0], replacement[1]);
            }
        }
        return key;
    }

    function reverseValue(key, value) {
        if (key === 'incline' && numeric.test(value)) {
            return value.replace(numeric, function(_, sign) { return sign === '-' ? '' : '-'; });
        } else if (key === 'incline' || key === 'direction') {
            return {up: 'down', down: 'up'}[value] || value;
        } else {
            return {left: 'right', right: 'left'}[value] || value;
        }
    }

    return function(graph) {
        var way = graph.entity(wayId),
            nodes = way.nodes.slice().reverse(),
            tags = {}, key, role;

        for (key in way.tags) {
            tags[reverseKey(key)] = reverseValue(key, way.tags[key]);
        }

        graph.parentRelations(way).forEach(function(relation) {
            relation.members.forEach(function(member, index) {
                if (member.id === way.id && (role = roleReversals[member.role])) {
                    relation = relation.updateMember({role: role}, index);
                    graph = graph.replace(relation);
                }
            });
        });

        return graph.replace(way.update({nodes: nodes, tags: tags}));
    };
};
iD.util = {};

iD.util.tagText = function(entity) {
    return d3.entries(entity.tags).map(function(e) {
        return e.key + '=' + e.value;
    }).join(', ');
};

iD.util.entitySelector = function(ids) {
    return ids.length ? '.' + ids.join(',.') : 'nothing';
};

iD.util.entityOrMemberSelector = function(ids, graph) {
    var s = iD.util.entitySelector(ids);

    ids.forEach(function(id) {
        var entity = graph.hasEntity(id);
        if (entity && entity.type === 'relation') {
            entity.members.forEach(function(member) {
                s += ',.' + member.id;
            });
        }
    });

    return s;
};

iD.util.displayName = function(entity) {
    var localeName = 'name:' + iD.detect().locale.toLowerCase().split('-')[0];
    return entity.tags[localeName] || entity.tags.name || entity.tags.ref;
};

iD.util.displayType = function(id) {
    return {
        n: t('inspector.node'),
        w: t('inspector.way'),
        r: t('inspector.relation')
    }[id.charAt(0)];
};

iD.util.stringQs = function(str) {
    return str.split('&').reduce(function(obj, pair){
        var parts = pair.split('=');
        if (parts.length === 2) {
            obj[parts[0]] = (null === parts[1]) ? '' : decodeURIComponent(parts[1]);
        }
        return obj;
    }, {});
};

iD.util.qsString = function(obj, noencode) {
    function softEncode(s) {
      // encode everything except special characters used in certain hash parameters:
      // "/" in map states, ":", ",", {" and "}" in background
      return encodeURIComponent(s).replace(/(%2F|%3A|%2C|%7B|%7D)/g, decodeURIComponent);
    }
    return Object.keys(obj).sort().map(function(key) {
        return encodeURIComponent(key) + '=' + (
            noencode ? softEncode(obj[key]) : encodeURIComponent(obj[key]));
    }).join('&');
};

iD.util.prefixDOMProperty = function(property) {
    var prefixes = ['webkit', 'ms', 'moz', 'o'],
        i = -1,
        n = prefixes.length,
        s = document.body;

    if (property in s)
        return property;

    property = property.substr(0, 1).toUpperCase() + property.substr(1);

    while (++i < n)
        if (prefixes[i] + property in s)
            return prefixes[i] + property;

    return false;
};

iD.util.prefixCSSProperty = function(property) {
    var prefixes = ['webkit', 'ms', 'Moz', 'O'],
        i = -1,
        n = prefixes.length,
        s = document.body.style;

    if (property.toLowerCase() in s)
        return property.toLowerCase();

    while (++i < n)
        if (prefixes[i] + property in s)
            return '-' + prefixes[i].toLowerCase() + property.replace(/([A-Z])/g, '-$1').toLowerCase();

    return false;
};


iD.util.setTransform = function(el, x, y, scale) {
    var prop = iD.util.transformProperty = iD.util.transformProperty || iD.util.prefixCSSProperty('Transform'),
        translate = iD.detect().opera ?
            'translate('   + x + 'px,' + y + 'px)' :
            'translate3d(' + x + 'px,' + y + 'px,0)';
    return el.style(prop, translate + (scale ? ' scale(' + scale + ')' : ''));
};

iD.util.getStyle = function(selector) {
    for (var i = 0; i < document.styleSheets.length; i++) {
        var rules = document.styleSheets[i].rules || document.styleSheets[i].cssRules || [];
        for (var k = 0; k < rules.length; k++) {
            var selectorText = rules[k].selectorText && rules[k].selectorText.split(', ');
            if (_.contains(selectorText, selector)) {
                return rules[k];
            }
        }
    }
};

iD.util.editDistance = function(a, b) {
    if (a.length === 0) return b.length;
    if (b.length === 0) return a.length;
    var matrix = [];
    for (var i = 0; i <= b.length; i++) { matrix[i] = [i]; }
    for (var j = 0; j <= a.length; j++) { matrix[0][j] = j; }
    for (i = 1; i <= b.length; i++) {
        for (j = 1; j <= a.length; j++) {
            if (b.charAt(i-1) === a.charAt(j-1)) {
                matrix[i][j] = matrix[i-1][j-1];
            } else {
                matrix[i][j] = Math.min(matrix[i-1][j-1] + 1, // substitution
                    Math.min(matrix[i][j-1] + 1, // insertion
                    matrix[i-1][j] + 1)); // deletion
            }
        }
    }
    return matrix[b.length][a.length];
};

// a d3.mouse-alike which
// 1. Only works on HTML elements, not SVG
// 2. Does not cause style recalculation
iD.util.fastMouse = function(container) {
    var rect = container.getBoundingClientRect(),
        rectLeft = rect.left,
        rectTop = rect.top,
        clientLeft = +container.clientLeft,
        clientTop = +container.clientTop;
    return function(e) {
        return [
            e.clientX - rectLeft - clientLeft,
            e.clientY - rectTop - clientTop];
    };
};

/* eslint-disable no-proto */
iD.util.getPrototypeOf = Object.getPrototypeOf || function(obj) { return obj.__proto__; };
/* eslint-enable no-proto */

iD.util.asyncMap = function(inputs, func, callback) {
    var remaining = inputs.length,
        results = [],
        errors = [];

    inputs.forEach(function(d, i) {
        func(d, function done(err, data) {
            errors[i] = err;
            results[i] = data;
            remaining--;
            if (!remaining) callback(errors, results);
        });
    });
};

// wraps an index to an interval [0..length-1]
iD.util.wrap = function(index, length) {
    if (index < 0)
        index += Math.ceil(-index/length)*length;
    return index % length;
};
iD.presets = function() {

    // an iD.presets.Collection with methods for
    // loading new data and returning defaults

    var all = iD.presets.Collection([]),
        defaults = { area: all, line: all, point: all, vertex: all, relation: all },
        fields = {},
        universal = [],
        recent = iD.presets.Collection([]);

    // Index of presets by (geometry, tag key).
    var index = {
        point: {},
        vertex: {},
        line: {},
        area: {},
        relation: {}
    };

    all.match = function(entity, resolver) {
        var geometry = entity.geometry(resolver),
            geometryMatches = index[geometry],
            best = -1,
            match;

        for (var k in entity.tags) {
            var keyMatches = geometryMatches[k];
            if (!keyMatches) continue;

            for (var i = 0; i < keyMatches.length; i++) {
                var score = keyMatches[i].matchScore(entity);
                if (score > best) {
                    best = score;
                    match = keyMatches[i];
                }
            }
        }

        return match || all.item(geometry);
    };

    // Because of the open nature of tagging, iD will never have a complete
    // list of tags used in OSM, so we want it to have logic like "assume
    // that a closed way with an amenity tag is an area, unless the amenity
    // is one of these specific types". This function computes a structure
    // that allows testing of such conditions, based on the presets designated
    // as as supporting (or not supporting) the area geometry.
    //
    // The returned object L is a whitelist/blacklist of tags. A closed way
    // with a tag (k, v) is considered to be an area if `k in L && !(v in L[k])`
    // (see `iD.Way#isArea()`). In other words, the keys of L form the whitelist,
    // and the subkeys form the blacklist.
    all.areaKeys = function() {
        var areaKeys = {},
            ignore = ['barrier', 'highway', 'footway', 'railway', 'type'],
            presets = _.reject(all.collection, 'suggestion');

        // whitelist
        presets.forEach(function(d) {
            for (var key in d.tags) break;
            if (!key) return;
            if (ignore.indexOf(key) !== -1) return;

            if (d.geometry.indexOf('area') !== -1) {
                areaKeys[key] = areaKeys[key] || {};
            }
        });

        // blacklist
        presets.forEach(function(d) {
            for (var key in d.tags) break;
            if (!key) return;
            if (ignore.indexOf(key) !== -1) return;

            var value = d.tags[key];
            if (d.geometry.indexOf('area') === -1 && key in areaKeys && value !== '*') {
                areaKeys[key][value] = true;
            }
        });

        return areaKeys;
    };

    all.load = function(d) {

        if (d.fields) {
            _.forEach(d.fields, function(d, id) {
                fields[id] = iD.presets.Field(id, d);
                if (d.universal) universal.push(fields[id]);
            });
        }

        if (d.presets) {
            _.forEach(d.presets, function(d, id) {
                all.collection.push(iD.presets.Preset(id, d, fields));
            });
        }

        if (d.categories) {
            _.forEach(d.categories, function(d, id) {
                all.collection.push(iD.presets.Category(id, d, all));
            });
        }

        if (d.defaults) {
            var getItem = _.bind(all.item, all);
            defaults = {
                area: iD.presets.Collection(d.defaults.area.map(getItem)),
                line: iD.presets.Collection(d.defaults.line.map(getItem)),
                point: iD.presets.Collection(d.defaults.point.map(getItem)),
                vertex: iD.presets.Collection(d.defaults.vertex.map(getItem)),
                relation: iD.presets.Collection(d.defaults.relation.map(getItem))
            };
        }

        for (var i = 0; i < all.collection.length; i++) {
            var preset = all.collection[i],
                geometry = preset.geometry;

            for (var j = 0; j < geometry.length; j++) {
                var g = index[geometry[j]];
                for (var k in preset.tags) {
                    (g[k] = g[k] || []).push(preset);
                }
            }
        }

        return all;
    };

    all.field = function(id) {
        return fields[id];
    };

    all.universal = function() {
        return universal;
    };

    all.defaults = function(geometry, n) {
        var rec = recent.matchGeometry(geometry).collection.slice(0, 4),
            def = _.uniq(rec.concat(defaults[geometry].collection)).slice(0, n - 1);
        return iD.presets.Collection(_.unique(rec.concat(def).concat(all.item(geometry))));
    };

    all.choose = function(preset) {
        if (!preset.isFallback()) {
            recent = iD.presets.Collection(_.unique([preset].concat(recent.collection)));
        }
        return all;
    };

    return all;
};
iD.presets.Category = function(id, category, all) {
    category = _.clone(category);

    category.id = id;

    category.members = iD.presets.Collection(category.members.map(function(id) {
        return all.item(id);
    }));

    category.matchGeometry = function(geometry) {
        return category.geometry.indexOf(geometry) >= 0;
    };

    category.matchScore = function() { return -1; };

    category.name = function() {
        return t('presets.categories.' + id + '.name', {'default': id});
    };

    category.terms = function() {
        return [];
    };

    return category;
};
iD.presets.Collection = function(collection) {

    var maxSearchResults = 50,
        maxSuggestionResults = 10;

    var presets = {

        collection: collection,

        item: function(id) {
            return _.find(collection, function(d) {
                return d.id === id;
            });
        },

        matchGeometry: function(geometry) {
            return iD.presets.Collection(collection.filter(function(d) {
                return d.matchGeometry(geometry);
            }));
        },

        search: function(value, geometry) {
            if (!value) return this;

            value = value.toLowerCase();

            var searchable = _.filter(collection, function(a) {
                    return a.searchable !== false && a.suggestion !== true;
                }),
                suggestions = _.filter(collection, function(a) {
                    return a.suggestion === true;
                });

            // matches value to preset.name
            var leading_name = _.filter(searchable, function(a) {
                    return leading(a.name().toLowerCase());
                }).sort(function(a, b) {
                    var i = a.name().toLowerCase().indexOf(value) - b.name().toLowerCase().indexOf(value);
                    if (i === 0) return a.name().length - b.name().length;
                    else return i;
                });

            // matches value to preset.terms values
            var leading_terms = _.filter(searchable, function(a) {
                    return _.any(a.terms() || [], leading);
                });

            // matches value to preset.tags values
            var leading_tag_values = _.filter(searchable, function(a) {
                    return _.any(_.without(_.values(a.tags || {}), '*'), leading);
                });

            function leading(a) {
                var index = a.indexOf(value);
                return index === 0 || a[index - 1] === ' ';
            }

            // finds close matches to value in preset.name
            var levenstein_name = searchable.map(function(a) {
                    return {
                        preset: a,
                        dist: iD.util.editDistance(value, a.name().toLowerCase())
                    };
                }).filter(function(a) {
                    return a.dist + Math.min(value.length - a.preset.name().length, 0) < 3;
                }).sort(function(a, b) {
                    return a.dist - b.dist;
                }).map(function(a) {
                    return a.preset;
                });

            // finds close matches to value in preset.terms
            var leventstein_terms = _.filter(searchable, function(a) {
                    return _.any(a.terms() || [], function(b) {
                        return iD.util.editDistance(value, b) + Math.min(value.length - b.length, 0) < 3;
                    });
                });

            function suggestionName(name) {
                var nameArray = name.split(' - ');
                if (nameArray.length > 1) {
                    name = nameArray.slice(0, nameArray.length-1).join(' - ');
                }
                return name.toLowerCase();
            }

            var leading_suggestions = _.filter(suggestions, function(a) {
                    return leading(suggestionName(a.name()));
                }).sort(function(a, b) {
                    a = suggestionName(a.name());
                    b = suggestionName(b.name());
                    var i = a.indexOf(value) - b.indexOf(value);
                    if (i === 0) return a.length - b.length;
                    else return i;
                });

            var leven_suggestions = suggestions.map(function(a) {
                    return {
                        preset: a,
                        dist: iD.util.editDistance(value, suggestionName(a.name()))
                    };
                }).filter(function(a) {
                    return a.dist + Math.min(value.length - suggestionName(a.preset.name()).length, 0) < 1;
                }).sort(function(a, b) {
                    return a.dist - b.dist;
                }).map(function(a) {
                    return a.preset;
                });

            var other = presets.item(geometry);

            var results = leading_name.concat(
                            leading_terms,
                            leading_tag_values,
                            leading_suggestions.slice(0, maxSuggestionResults+5),
                            levenstein_name,
                            leventstein_terms,
                            leven_suggestions.slice(0, maxSuggestionResults)
                        ).slice(0, maxSearchResults-1);

            return iD.presets.Collection(_.unique(
                    results.concat(other)
                ));
        }
    };

    return presets;
};
iD.presets.Field = function(id, field) {
    field = _.clone(field);

    field.id = id;

    field.matchGeometry = function(geometry) {
        return !field.geometry || field.geometry === geometry;
    };

    field.t = function(scope, options) {
        return t('presets.fields.' + id + '.' + scope, options);
    };

    field.label = function() {
        return field.t('label', {'default': id});
    };

    var placeholder = field.placeholder;
    field.placeholder = function() {
        return field.t('placeholder', {'default': placeholder});
    };

    return field;
};
iD.presets.Preset = function(id, preset, fields) {
    preset = _.clone(preset);

    preset.id = id;
    preset.fields = (preset.fields || []).map(getFields);
    preset.geometry = (preset.geometry || []);

    function getFields(f) {
        return fields[f];
    }

    preset.matchGeometry = function(geometry) {
        return preset.geometry.indexOf(geometry) >= 0;
    };

    var matchScore = preset.matchScore || 1;
    preset.matchScore = function(entity) {
        var tags = preset.tags,
            score = 0;

        for (var t in tags) {
            if (entity.tags[t] === tags[t]) {
                score += matchScore;
            } else if (tags[t] === '*' && t in entity.tags) {
                score += matchScore / 2;
            } else {
                return -1;
            }
        }

        return score;
    };

    preset.t = function(scope, options) {
        return t('presets.presets.' + id + '.' + scope, options);
    };

    var name = preset.name;
    preset.name = function() {
        if (preset.suggestion) {
            id = id.split('/');
            id = id[0] + '/' + id[1];
            return name + ' - ' + t('presets.presets.' + id + '.name');
        }
        return preset.t('name', {'default': name});
    };

    preset.terms = function() {
        return preset.t('terms', {'default': ''}).toLowerCase().split(/\s*,+\s*/);
    };

    preset.isFallback = function() {
        return Object.keys(preset.tags).length === 0;
    };

    preset.reference = function(geometry) {
        var key = Object.keys(preset.tags)[0],
            value = preset.tags[key];

        if (geometry === 'relation' && key === 'type') {
            return { rtype: value };
        } else if (value === '*') {
            return { key: key };
        } else {
            return { key: key, value: value };
        }
    };

    var removeTags = preset.removeTags || preset.tags;
    preset.removeTags = function(tags, geometry) {
        tags = _.omit(tags, _.keys(removeTags));

        for (var f in preset.fields) {
            var field = preset.fields[f];
            if (field.matchGeometry(geometry) && field.default === tags[field.key]) {
                delete tags[field.key];
            }
        }

        delete tags.area;
        return tags;
    };

    var applyTags = preset.addTags || preset.tags;
    preset.applyTags = function(tags, geometry) {
        var k;

        tags = _.clone(tags);

        for (k in applyTags) {
            if (applyTags[k] === '*') {
                tags[k] = 'yes';
            } else {
                tags[k] = applyTags[k];
            }
        }

        // Add area=yes if necessary.
        // This is necessary if the geometry is already an area (e.g. user drew an area) AND any of:
        // 1. chosen preset could be either an area or a line (`barrier=city_wall`)
        // 2. chosen preset doesn't have a key in areaKeys (`railway=station`)
        if (geometry === 'area') {
            var needsAreaTag = true;
            if (preset.geometry.indexOf('line') === -1) {
                for (k in applyTags) {
                    if (k in iD.areaKeys) {
                        needsAreaTag = false;
                        break;
                    }
                }
            }
            if (needsAreaTag) {
                tags.area = 'yes';
            }
        }

        for (var f in preset.fields) {
            var field = preset.fields[f];
            if (field.matchGeometry(geometry) && field.key && !tags[field.key] && field.default) {
                tags[field.key] = field.default;
            }
        }

        return tags;
    };

    return preset;
};
/*
    iD.Difference represents the difference between two graphs.
    It knows how to calculate the set of entities that were
    created, modified, or deleted, and also contains the logic
    for recursively extending a difference to the complete set
    of entities that will require a redraw, taking into account
    child and parent relationships.
 */
iD.Difference = function(base, head) {
    var changes = {}, length = 0;

    function changed(h, b) {
        return h !== b && !_.isEqual(_.omit(h, 'v'), _.omit(b, 'v'));
    }

    _.each(head.entities, function(h, id) {
        var b = base.entities[id];
        if (changed(h, b)) {
            changes[id] = {base: b, head: h};
            length++;
        }
    });

    _.each(base.entities, function(b, id) {
        var h = head.entities[id];
        if (!changes[id] && changed(h, b)) {
            changes[id] = {base: b, head: h};
            length++;
        }
    });

    function addParents(parents, result) {
        for (var i = 0; i < parents.length; i++) {
            var parent = parents[i];

            if (parent.id in result)
                continue;

            result[parent.id] = parent;
            addParents(head.parentRelations(parent), result);
        }
    }

    var difference = {};

    difference.length = function() {
        return length;
    };

    difference.changes = function() {
        return changes;
    };

    difference.extantIDs = function() {
        var result = [];
        _.each(changes, function(change, id) {
            if (change.head) result.push(id);
        });
        return result;
    };

    difference.modified = function() {
        var result = [];
        _.each(changes, function(change) {
            if (change.base && change.head) result.push(change.head);
        });
        return result;
    };

    difference.created = function() {
        var result = [];
        _.each(changes, function(change) {
            if (!change.base && change.head) result.push(change.head);
        });
        return result;
    };

    difference.deleted = function() {
        var result = [];
        _.each(changes, function(change) {
            if (change.base && !change.head) result.push(change.base);
        });
        return result;
    };

    difference.summary = function() {
        var relevant = {};

        function addEntity(entity, graph, changeType) {
            relevant[entity.id] = {
                entity: entity,
                graph: graph,
                changeType: changeType
            };
        }

        function addParents(entity) {
            var parents = head.parentWays(entity);
            for (var j = parents.length - 1; j >= 0; j--) {
                var parent = parents[j];
                if (!(parent.id in relevant)) addEntity(parent, head, 'modified');
            }
        }

        _.each(changes, function(change) {
            if (change.head && change.head.geometry(head) !== 'vertex') {
                addEntity(change.head, head, change.base ? 'modified' : 'created');

            } else if (change.base && change.base.geometry(base) !== 'vertex') {
                addEntity(change.base, base, 'deleted');

            } else if (change.base && change.head) { // modified vertex
                var moved    = !_.isEqual(change.base.loc,  change.head.loc),
                    retagged = !_.isEqual(change.base.tags, change.head.tags);

                if (moved) {
                    addParents(change.head);
                }

                if (retagged || (moved && change.head.hasInterestingTags())) {
                    addEntity(change.head, head, 'modified');
                }

            } else if (change.head && change.head.hasInterestingTags()) { // created vertex
                addEntity(change.head, head, 'created');

            } else if (change.base && change.base.hasInterestingTags()) { // deleted vertex
                addEntity(change.base, base, 'deleted');
            }
        });

        return d3.values(relevant);
    };

    difference.complete = function(extent) {
        var result = {}, id, change;

        for (id in changes) {
            change = changes[id];

            var h = change.head,
                b = change.base,
                entity = h || b;

            if (extent &&
                (!h || !h.intersects(extent, head)) &&
                (!b || !b.intersects(extent, base)))
                continue;

            result[id] = h;

            if (entity.type === 'way') {
                var nh = h ? h.nodes : [],
                    nb = b ? b.nodes : [],
                    diff, i;

                diff = _.difference(nh, nb);
                for (i = 0; i < diff.length; i++) {
                    result[diff[i]] = head.hasEntity(diff[i]);
                }

                diff = _.difference(nb, nh);
                for (i = 0; i < diff.length; i++) {
                    result[diff[i]] = head.hasEntity(diff[i]);
                }
            }

            addParents(head.parentWays(entity), result);
            addParents(head.parentRelations(entity), result);
        }

        return result;
    };

    return difference;
};
iD.Tree = function(head) {
    var rtree = rbush(),
        rectangles = {};

    function extentRectangle(extent) {
        return [
            extent[0][0],
            extent[0][1],
            extent[1][0],
            extent[1][1]
        ];
    }

    function entityRectangle(entity) {
        var rect = extentRectangle(entity.extent(head));
        rect.id = entity.id;
        rectangles[entity.id] = rect;
        return rect;
    }

    function updateParents(entity, insertions, memo) {
        head.parentWays(entity).forEach(function(parent) {
            if (rectangles[parent.id]) {
                rtree.remove(rectangles[parent.id]);
                insertions[parent.id] = parent;
            }
        });

        head.parentRelations(entity).forEach(function(parent) {
            if (memo[entity.id]) return;
            memo[entity.id] = true;
            if (rectangles[parent.id]) {
                rtree.remove(rectangles[parent.id]);
                insertions[parent.id] = parent;
            }
            updateParents(parent, insertions, memo);
        });
    }

    var tree = {};

    tree.rebase = function(entities, force) {
        var insertions = {};

        for (var i = 0; i < entities.length; i++) {
            var entity = entities[i];

            if (!entity.visible)
                continue;

            if (head.entities.hasOwnProperty(entity.id) || rectangles[entity.id]) {
                if (!force) {
                    continue;
                } else if (rectangles[entity.id]) {
                    rtree.remove(rectangles[entity.id]);
                }
            }

            insertions[entity.id] = entity;
            updateParents(entity, insertions, {});
        }

        rtree.load(_.map(insertions, entityRectangle));

        return tree;
    };

    tree.intersects = function(extent, graph) {
        if (graph !== head) {
            var diff = iD.Difference(head, graph),
                insertions = {};

            head = graph;

            diff.deleted().forEach(function(entity) {
                rtree.remove(rectangles[entity.id]);
                delete rectangles[entity.id];
            });

            diff.modified().forEach(function(entity) {
                rtree.remove(rectangles[entity.id]);
                insertions[entity.id] = entity;
                updateParents(entity, insertions, {});
            });

            diff.created().forEach(function(entity) {
                insertions[entity.id] = entity;
            });

            rtree.load(_.map(insertions, entityRectangle));
        }

        return rtree.search(extentRectangle(extent)).map(function(rect) {
            return head.entity(rect.id);
        });
    };

    return tree;
};
iD.Entity = function(attrs) {
    // For prototypal inheritance.
    if (this instanceof iD.Entity) return;

    // Create the appropriate subtype.
    if (attrs && attrs.type) {
        return iD.Entity[attrs.type].apply(this, arguments);
    } else if (attrs && attrs.id) {
        return iD.Entity[iD.Entity.id.type(attrs.id)].apply(this, arguments);
    }

    // Initialize a generic Entity (used only in tests).
    return (new iD.Entity()).initialize(arguments);
};

iD.Entity.id = function(type) {
    return iD.Entity.id.fromOSM(type, iD.Entity.id.next[type]--);
};

iD.Entity.id.next = {node: -1, way: -1, relation: -1};

iD.Entity.id.fromOSM = function(type, id) {
    return type[0] + id;
};

iD.Entity.id.toOSM = function(id) {
    return id.slice(1);
};

iD.Entity.id.type = function(id) {
    return {'n': 'node', 'w': 'way', 'r': 'relation'}[id[0]];
};

// A function suitable for use as the second argument to d3.selection#data().
iD.Entity.key = function(entity) {
    return entity.id + 'v' + (entity.v || 0);
};

iD.Entity.prototype = {
    tags: {},

    initialize: function(sources) {
        for (var i = 0; i < sources.length; ++i) {
            var source = sources[i];
            for (var prop in source) {
                if (Object.prototype.hasOwnProperty.call(source, prop)) {
                    if (source[prop] === undefined) {
                        delete this[prop];
                    } else {
                        this[prop] = source[prop];
                    }
                }
            }
        }

        if (!this.id && this.type) {
            this.id = iD.Entity.id(this.type);
        }
        if (!this.hasOwnProperty('visible')) {
            this.visible = true;
        }

        if (iD.debug) {
            Object.freeze(this);
            Object.freeze(this.tags);

            if (this.loc) Object.freeze(this.loc);
            if (this.nodes) Object.freeze(this.nodes);
            if (this.members) Object.freeze(this.members);
        }

        return this;
    },

    copy: function() {
        // Returns an array so that we can support deep copying ways and relations.
        // The first array element will contain this.copy, followed by any descendants.
        return [iD.Entity(this, {id: undefined, user: undefined, version: undefined})];
    },

    osmId: function() {
        return iD.Entity.id.toOSM(this.id);
    },

    isNew: function() {
        return this.osmId() < 0;
    },

    update: function(attrs) {
        return iD.Entity(this, attrs, {v: 1 + (this.v || 0)});
    },

    mergeTags: function(tags) {
        var merged = _.clone(this.tags), changed = false;
        for (var k in tags) {
            var t1 = merged[k],
                t2 = tags[k];
            if (!t1) {
                changed = true;
                merged[k] = t2;
            } else if (t1 !== t2) {
                changed = true;
                merged[k] = _.union(t1.split(/;\s*/), t2.split(/;\s*/)).join(';');
            }
        }
        return changed ? this.update({tags: merged}) : this;
    },

    intersects: function(extent, resolver) {
        return this.extent(resolver).intersects(extent);
    },

    isUsed: function(resolver) {
        return _.without(Object.keys(this.tags), 'area').length > 0 ||
            resolver.parentRelations(this).length > 0;
    },

    hasInterestingTags: function() {
        return _.keys(this.tags).some(function(key) {
            return key !== 'attribution' &&
                key !== 'created_by' &&
                key !== 'source' &&
                key !== 'odbl' &&
                key.indexOf('tiger:') !== 0;
        });
    },

    isHighwayIntersection: function() {
        return false;
    },

    deprecatedTags: function() {
        var tags = _.pairs(this.tags);
        var deprecated = {};

        iD.data.deprecated.forEach(function(d) {
            var match = _.pairs(d.old)[0];
            tags.forEach(function(t) {
                if (t[0] === match[0] &&
                    (t[1] === match[1] || match[1] === '*')) {
                    deprecated[t[0]] = t[1];
                }
            });
        });

        return deprecated;
    }
};
iD.Node = iD.Entity.node = function iD_Node() {
    if (!(this instanceof iD_Node)) {
        return (new iD_Node()).initialize(arguments);
    } else if (arguments.length) {
        this.initialize(arguments);
    }
};

iD.Node.prototype = Object.create(iD.Entity.prototype);

_.extend(iD.Node.prototype, {
    type: 'node',

    extent: function() {
        return new iD.geo.Extent(this.loc);
    },

    geometry: function(graph) {
        return graph.transient(this, 'geometry', function() {
            return graph.isPoi(this) ? 'point' : 'vertex';
        });
    },

    move: function(loc) {
        return this.update({loc: loc});
    },

    isIntersection: function(resolver) {
        return resolver.transient(this, 'isIntersection', function() {
            return resolver.parentWays(this).filter(function(parent) {
                return (parent.tags.highway ||
                    parent.tags.waterway ||
                    parent.tags.railway ||
                    parent.tags.aeroway) &&
                    parent.geometry(resolver) === 'line';
            }).length > 1;
        });
    },

    isHighwayIntersection: function(resolver) {
        return resolver.transient(this, 'isHighwayIntersection', function() {
            return resolver.parentWays(this).filter(function(parent) {
                return parent.tags.highway && parent.geometry(resolver) === 'line';
            }).length > 1;
        });
    },

    asJXON: function(changeset_id) {
        var r = {
            node: {
                '@id': this.osmId(),
                '@lon': this.loc[0],
                '@lat': this.loc[1],
                '@version': (this.version || 0),
                tag: _.map(this.tags, function(v, k) {
                    return { keyAttributes: { k: k, v: v } };
                })
            }
        };
        if (changeset_id) r.node['@changeset'] = changeset_id;
        return r;
    },

    asGeoJSON: function() {
        return {
            type: 'Point',
            coordinates: this.loc
        };
    }
});
iD.Way = iD.Entity.way = function iD_Way() {
    if (!(this instanceof iD_Way)) {
        return (new iD_Way()).initialize(arguments);
    } else if (arguments.length) {
        this.initialize(arguments);
    }
};

iD.Way.prototype = Object.create(iD.Entity.prototype);

_.extend(iD.Way.prototype, {
    type: 'way',
    nodes: [],

    copy: function(deep, resolver) {
        var copy = iD.Entity.prototype.copy.call(this);

        if (!deep || !resolver) {
            return copy;
        }

        var nodes = [],
            replacements = {},
            i, oldid, newid, child;

        for (i = 0; i < this.nodes.length; i++) {
            oldid = this.nodes[i];
            newid = replacements[oldid];
            if (!newid) {
                child = resolver.entity(oldid).copy();
                newid = replacements[oldid] = child[0].id;
                copy = copy.concat(child);
            }
            nodes.push(newid);
        }

        copy[0] = copy[0].update({nodes: nodes});
        return copy;
    },

    extent: function(resolver) {
        return resolver.transient(this, 'extent', function() {
            var extent = iD.geo.Extent();
            for (var i = 0; i < this.nodes.length; i++) {
                var node = resolver.hasEntity(this.nodes[i]);
                if (node) {
                    extent._extend(node.extent());
                }
            }
            return extent;
        });
    },

    first: function() {
        return this.nodes[0];
    },

    last: function() {
        return this.nodes[this.nodes.length - 1];
    },

    contains: function(node) {
        return this.nodes.indexOf(node) >= 0;
    },

    affix: function(node) {
        if (this.nodes[0] === node) return 'prefix';
        if (this.nodes[this.nodes.length - 1] === node) return 'suffix';
    },

    layer: function() {
        // explicit layer tag, clamp between -10, 10..
        if (this.tags.layer !== undefined) {
            return Math.max(-10, Math.min(+(this.tags.layer), 10));
        }

        // implied layer tag..
        if (this.tags.location === 'overground') return 1;
        if (this.tags.location === 'underground') return -1;
        if (this.tags.location === 'underwater') return -10;

        if (this.tags.power === 'line') return 10;
        if (this.tags.power === 'minor_line') return 10;
        if (this.tags.aerialway) return 10;
        if (this.tags.bridge) return 1;
        if (this.tags.cutting) return -1;
        if (this.tags.tunnel) return -1;
        if (this.tags.waterway) return -1;
        if (this.tags.man_made === 'pipeline') return -10;
        if (this.tags.boundary) return -10;
        return 0;
    },

    isOneWay: function() {
        // explicit oneway tag..
        if (['yes', '1', '-1'].indexOf(this.tags.oneway) !== -1) { return true; }
        if (['no', '0'].indexOf(this.tags.oneway) !== -1) { return false; }

        // implied oneway tag..
        for (var key in this.tags) {
            if (key in iD.oneWayTags && (this.tags[key] in iD.oneWayTags[key]))
                return true;
        }
        return false;
    },

    isClosed: function() {
        return this.nodes.length > 0 && this.first() === this.last();
    },

    isConvex: function(resolver) {
        if (!this.isClosed() || this.isDegenerate()) return null;

        var nodes = _.uniq(resolver.childNodes(this)),
            coords = _.pluck(nodes, 'loc'),
            curr = 0, prev = 0;

        for (var i = 0; i < coords.length; i++) {
            var o = coords[(i+1) % coords.length],
                a = coords[i],
                b = coords[(i+2) % coords.length],
                res = iD.geo.cross(o, a, b);

            curr = (res > 0) ? 1 : (res < 0) ? -1 : 0;
            if (curr === 0) {
                continue;
            } else if (prev && curr !== prev) {
                return false;
            }
            prev = curr;
        }
        return true;
    },

    isArea: function() {
        if (this.tags.area === 'yes')
            return true;
        if (!this.isClosed() || this.tags.area === 'no')
            return false;
        for (var key in this.tags)
            if (key in iD.areaKeys && !(this.tags[key] in iD.areaKeys[key]))
                return true;
        return false;
    },

    isDegenerate: function() {
        return _.uniq(this.nodes).length < (this.isArea() ? 3 : 2);
    },

    areAdjacent: function(n1, n2) {
        for (var i = 0; i < this.nodes.length; i++) {
            if (this.nodes[i] === n1) {
                if (this.nodes[i - 1] === n2) return true;
                if (this.nodes[i + 1] === n2) return true;
            }
        }
        return false;
    },

    geometry: function(graph) {
        return graph.transient(this, 'geometry', function() {
            return this.isArea() ? 'area' : 'line';
        });
    },

    addNode: function(id, index) {
        var nodes = this.nodes.slice();
        nodes.splice(index === undefined ? nodes.length : index, 0, id);
        return this.update({nodes: nodes});
    },

    updateNode: function(id, index) {
        var nodes = this.nodes.slice();
        nodes.splice(index, 1, id);
        return this.update({nodes: nodes});
    },

    replaceNode: function(needle, replacement) {
        if (this.nodes.indexOf(needle) < 0)
            return this;

        var nodes = this.nodes.slice();
        for (var i = 0; i < nodes.length; i++) {
            if (nodes[i] === needle) {
                nodes[i] = replacement;
            }
        }
        return this.update({nodes: nodes});
    },

    removeNode: function(id) {
        var nodes = [];

        for (var i = 0; i < this.nodes.length; i++) {
            var node = this.nodes[i];
            if (node !== id && nodes[nodes.length - 1] !== node) {
                nodes.push(node);
            }
        }

        // Preserve circularity
        if (this.nodes.length > 1 && this.first() === id && this.last() === id && nodes[nodes.length - 1] !== nodes[0]) {
            nodes.push(nodes[0]);
        }

        return this.update({nodes: nodes});
    },

    asJXON: function(changeset_id) {
        var r = {
            way: {
                '@id': this.osmId(),
                '@version': this.version || 0,
                nd: _.map(this.nodes, function(id) {
                    return { keyAttributes: { ref: iD.Entity.id.toOSM(id) } };
                }),
                tag: _.map(this.tags, function(v, k) {
                    return { keyAttributes: { k: k, v: v } };
                })
            }
        };
        if (changeset_id) r.way['@changeset'] = changeset_id;
        return r;
    },

    asGeoJSON: function(resolver) {
        return resolver.transient(this, 'GeoJSON', function() {
            var coordinates = _.pluck(resolver.childNodes(this), 'loc');
            if (this.isArea() && this.isClosed()) {
                return {
                    type: 'Polygon',
                    coordinates: [coordinates]
                };
            } else {
                return {
                    type: 'LineString',
                    coordinates: coordinates
                };
            }
        });
    },

    area: function(resolver) {
        return resolver.transient(this, 'area', function() {
            var nodes = resolver.childNodes(this);

            var json = {
                type: 'Polygon',
                coordinates: [_.pluck(nodes, 'loc')]
            };

            if (!this.isClosed() && nodes.length) {
                json.coordinates[0].push(nodes[0].loc);
            }

            var area = d3.geo.area(json);

            // Heuristic for detecting counterclockwise winding order. Assumes
            // that OpenStreetMap polygons are not hemisphere-spanning.
            if (area > 2 * Math.PI) {
                json.coordinates[0] = json.coordinates[0].reverse();
                area = d3.geo.area(json);
            }

            return isNaN(area) ? 0 : area;
        });
    }
});
iD.Relation = iD.Entity.relation = function iD_Relation() {
    if (!(this instanceof iD_Relation)) {
        return (new iD_Relation()).initialize(arguments);
    } else if (arguments.length) {
        this.initialize(arguments);
    }
};

iD.Relation.prototype = Object.create(iD.Entity.prototype);

iD.Relation.creationOrder = function(a, b) {
    var aId = parseInt(iD.Entity.id.toOSM(a.id), 10);
    var bId = parseInt(iD.Entity.id.toOSM(b.id), 10);

    if (aId < 0 || bId < 0) return aId - bId;
    return bId - aId;
};

_.extend(iD.Relation.prototype, {
    type: 'relation',
    members: [],

    copy: function(deep, resolver, replacements) {
        var copy = iD.Entity.prototype.copy.call(this);
        if (!deep || !resolver || !this.isComplete(resolver)) {
            return copy;
        }

        var members = [],
            i, oldmember, oldid, newid, children;

        replacements = replacements || {};
        replacements[this.id] = copy[0].id;

        for (i = 0; i < this.members.length; i++) {
            oldmember = this.members[i];
            oldid = oldmember.id;
            newid = replacements[oldid];
            if (!newid) {
                children = resolver.entity(oldid).copy(true, resolver, replacements);
                newid = replacements[oldid] = children[0].id;
                copy = copy.concat(children);
            }
            members.push({id: newid, type: oldmember.type, role: oldmember.role});
        }

        copy[0] = copy[0].update({members: members});
        return copy;
    },

    extent: function(resolver, memo) {
        return resolver.transient(this, 'extent', function() {
            if (memo && memo[this.id]) return iD.geo.Extent();
            memo = memo || {};
            memo[this.id] = true;

            var extent = iD.geo.Extent();
            for (var i = 0; i < this.members.length; i++) {
                var member = resolver.hasEntity(this.members[i].id);
                if (member) {
                    extent._extend(member.extent(resolver, memo));
                }
            }
            return extent;
        });
    },

    geometry: function(graph) {
        return graph.transient(this, 'geometry', function() {
            return this.isMultipolygon() ? 'area' : 'relation';
        });
    },

    isDegenerate: function() {
        return this.members.length === 0;
    },

    // Return an array of members, each extended with an 'index' property whose value
    // is the member index.
    indexedMembers: function() {
        var result = new Array(this.members.length);
        for (var i = 0; i < this.members.length; i++) {
            result[i] = _.extend({}, this.members[i], {index: i});
        }
        return result;
    },

    // Return the first member with the given role. A copy of the member object
    // is returned, extended with an 'index' property whose value is the member index.
    memberByRole: function(role) {
        for (var i = 0; i < this.members.length; i++) {
            if (this.members[i].role === role) {
                return _.extend({}, this.members[i], {index: i});
            }
        }
    },

    // Return the first member with the given id. A copy of the member object
    // is returned, extended with an 'index' property whose value is the member index.
    memberById: function(id) {
        for (var i = 0; i < this.members.length; i++) {
            if (this.members[i].id === id) {
                return _.extend({}, this.members[i], {index: i});
            }
        }
    },

    // Return the first member with the given id and role. A copy of the member object
    // is returned, extended with an 'index' property whose value is the member index.
    memberByIdAndRole: function(id, role) {
        for (var i = 0; i < this.members.length; i++) {
            if (this.members[i].id === id && this.members[i].role === role) {
                return _.extend({}, this.members[i], {index: i});
            }
        }
    },

    addMember: function(member, index) {
        var members = this.members.slice();
        members.splice(index === undefined ? members.length : index, 0, member);
        return this.update({members: members});
    },

    updateMember: function(member, index) {
        var members = this.members.slice();
        members.splice(index, 1, _.extend({}, members[index], member));
        return this.update({members: members});
    },

    removeMember: function(index) {
        var members = this.members.slice();
        members.splice(index, 1);
        return this.update({members: members});
    },

    removeMembersWithID: function(id) {
        var members = _.reject(this.members, function(m) { return m.id === id; });
        return this.update({members: members});
    },

    // Wherever a member appears with id `needle.id`, replace it with a member
    // with id `replacement.id`, type `replacement.type`, and the original role,
    // unless a member already exists with that id and role. Return an updated
    // relation.
    replaceMember: function(needle, replacement) {
        if (!this.memberById(needle.id))
            return this;

        var members = [];

        for (var i = 0; i < this.members.length; i++) {
            var member = this.members[i];
            if (member.id !== needle.id) {
                members.push(member);
            } else if (!this.memberByIdAndRole(replacement.id, member.role)) {
                members.push({id: replacement.id, type: replacement.type, role: member.role});
            }
        }

        return this.update({members: members});
    },

    asJXON: function(changeset_id) {
        var r = {
            relation: {
                '@id': this.osmId(),
                '@version': this.version || 0,
                member: _.map(this.members, function(member) {
                    return { keyAttributes: { type: member.type, role: member.role, ref: iD.Entity.id.toOSM(member.id) } };
                }),
                tag: _.map(this.tags, function(v, k) {
                    return { keyAttributes: { k: k, v: v } };
                })
            }
        };
        if (changeset_id) r.relation['@changeset'] = changeset_id;
        return r;
    },

    asGeoJSON: function(resolver) {
        return resolver.transient(this, 'GeoJSON', function () {
            if (this.isMultipolygon()) {
                return {
                    type: 'MultiPolygon',
                    coordinates: this.multipolygon(resolver)
                };
            } else {
                return {
                    type: 'FeatureCollection',
                    properties: this.tags,
                    features: this.members.map(function (member) {
                        return _.extend({role: member.role}, resolver.entity(member.id).asGeoJSON(resolver));
                    })
                };
            }
        });
    },

    area: function(resolver) {
        return resolver.transient(this, 'area', function() {
            return d3.geo.area(this.asGeoJSON(resolver));
        });
    },

    isMultipolygon: function() {
        return this.tags.type === 'multipolygon';
    },

    isComplete: function(resolver) {
        for (var i = 0; i < this.members.length; i++) {
            if (!resolver.hasEntity(this.members[i].id)) {
                return false;
            }
        }
        return true;
    },

    isRestriction: function() {
        return !!(this.tags.type && this.tags.type.match(/^restriction:?/));
    },

    // Returns an array [A0, ... An], each Ai being an array of node arrays [Nds0, ... Ndsm],
    // where Nds0 is an outer ring and subsequent Ndsi's (if any i > 0) being inner rings.
    //
    // This corresponds to the structure needed for rendering a multipolygon path using a
    // `evenodd` fill rule, as well as the structure of a GeoJSON MultiPolygon geometry.
    //
    // In the case of invalid geometries, this function will still return a result which
    // includes the nodes of all way members, but some Nds may be unclosed and some inner
    // rings not matched with the intended outer ring.
    //
    multipolygon: function(resolver) {
        var outers = this.members.filter(function(m) { return 'outer' === (m.role || 'outer'); }),
            inners = this.members.filter(function(m) { return 'inner' === m.role; });

        outers = iD.geo.joinWays(outers, resolver);
        inners = iD.geo.joinWays(inners, resolver);

        outers = outers.map(function(outer) { return _.pluck(outer.nodes, 'loc'); });
        inners = inners.map(function(inner) { return _.pluck(inner.nodes, 'loc'); });

        var result = outers.map(function(o) {
            // Heuristic for detecting counterclockwise winding order. Assumes
            // that OpenStreetMap polygons are not hemisphere-spanning.
            return [d3.geo.area({type: 'Polygon', coordinates: [o]}) > 2 * Math.PI ? o.reverse() : o];
        });

        function findOuter(inner) {
            var o, outer;

            for (o = 0; o < outers.length; o++) {
                outer = outers[o];
                if (iD.geo.polygonContainsPolygon(outer, inner))
                    return o;
            }

            for (o = 0; o < outers.length; o++) {
                outer = outers[o];
                if (iD.geo.polygonIntersectsPolygon(outer, inner))
                    return o;
            }
        }

        for (var i = 0; i < inners.length; i++) {
            var inner = inners[i];

            if (d3.geo.area({type: 'Polygon', coordinates: [inner]}) < 2 * Math.PI) {
                inner = inner.reverse();
            }

            var o = findOuter(inners[i]);
            if (o !== undefined)
                result[o].push(inners[i]);
            else
                result.push([inners[i]]); // Invalid geometry
        }

        return result;
    }
});
iD.Graph = function(other, mutable) {
    if (!(this instanceof iD.Graph)) return new iD.Graph(other, mutable);

    if (other instanceof iD.Graph) {
        var base = other.base();
        this.entities = _.assign(Object.create(base.entities), other.entities);
        this._parentWays = _.assign(Object.create(base.parentWays), other._parentWays);
        this._parentRels = _.assign(Object.create(base.parentRels), other._parentRels);

    } else {
        this.entities = Object.create({});
        this._parentWays = Object.create({});
        this._parentRels = Object.create({});
        this.rebase(other || [], [this]);
    }

    this.transients = {};
    this._childNodes = {};
    this.frozen = !mutable;
};

iD.Graph.prototype = {
    hasEntity: function(id) {
        return this.entities[id];
    },

    entity: function(id) {
        var entity = this.entities[id];
        if (!entity) {
            throw new Error('entity ' + id + ' not found');
        }
        return entity;
    },

    transient: function(entity, key, fn) {
        var id = entity.id,
            transients = this.transients[id] ||
            (this.transients[id] = {});

        if (transients[key] !== undefined) {
            return transients[key];
        }

        transients[key] = fn.call(entity);

        return transients[key];
    },

    parentWays: function(entity) {
        var parents = this._parentWays[entity.id],
            result = [];

        if (parents) {
            for (var i = 0; i < parents.length; i++) {
                result.push(this.entity(parents[i]));
            }
        }
        return result;
    },

    isPoi: function(entity) {
        var parentWays = this._parentWays[entity.id];
        return !parentWays || parentWays.length === 0;
    },

    isShared: function(entity) {
        var parentWays = this._parentWays[entity.id];
        return parentWays && parentWays.length > 1;
    },

    parentRelations: function(entity) {
        var parents = this._parentRels[entity.id],
            result = [];

        if (parents) {
            for (var i = 0; i < parents.length; i++) {
                result.push(this.entity(parents[i]));
            }
        }
        return result;
    },

    childNodes: function(entity) {
        if (this._childNodes[entity.id])
            return this._childNodes[entity.id];

        var nodes = [];
        if (entity.nodes) {
            for (var i = 0; i < entity.nodes.length; i++) {
                nodes[i] = this.entity(entity.nodes[i]);
            }
        }

        if (iD.debug) Object.freeze(nodes);

        this._childNodes[entity.id] = nodes;
        return this._childNodes[entity.id];
    },

    base: function() {
        return {
            'entities': iD.util.getPrototypeOf(this.entities),
            'parentWays': iD.util.getPrototypeOf(this._parentWays),
            'parentRels': iD.util.getPrototypeOf(this._parentRels)
        };
    },

    // Unlike other graph methods, rebase mutates in place. This is because it
    // is used only during the history operation that merges newly downloaded
    // data into each state. To external consumers, it should appear as if the
    // graph always contained the newly downloaded data.
    rebase: function(entities, stack, force) {
        var base = this.base(),
            i, j, k, id;

        for (i = 0; i < entities.length; i++) {
            var entity = entities[i];

            if (!entity.visible || (!force && base.entities[entity.id]))
                continue;

            // Merging data into the base graph
            base.entities[entity.id] = entity;
            this._updateCalculated(undefined, entity, base.parentWays, base.parentRels);

            // Restore provisionally-deleted nodes that are discovered to have an extant parent
            if (entity.type === 'way') {
                for (j = 0; j < entity.nodes.length; j++) {
                    id = entity.nodes[j];
                    for (k = 1; k < stack.length; k++) {
                        var ents = stack[k].entities;
                        if (ents.hasOwnProperty(id) && ents[id] === undefined) {
                            delete ents[id];
                        }
                    }
                }
            }
        }

        for (i = 0; i < stack.length; i++) {
            stack[i]._updateRebased();
        }
    },

    _updateRebased: function() {
        var base = this.base(),
            i, k, child, id, keys;

        keys = Object.keys(this._parentWays);
        for (i = 0; i < keys.length; i++) {
            child = keys[i];
            if (base.parentWays[child]) {
                for (k = 0; k < base.parentWays[child].length; k++) {
                    id = base.parentWays[child][k];
                    if (!this.entities.hasOwnProperty(id) && !_.contains(this._parentWays[child], id)) {
                        this._parentWays[child].push(id);
                    }
                }
            }
        }

        keys = Object.keys(this._parentRels);
        for (i = 0; i < keys.length; i++) {
            child = keys[i];
            if (base.parentRels[child]) {
                for (k = 0; k < base.parentRels[child].length; k++) {
                    id = base.parentRels[child][k];
                    if (!this.entities.hasOwnProperty(id) && !_.contains(this._parentRels[child], id)) {
                        this._parentRels[child].push(id);
                    }
                }
            }
        }

        this.transients = {};

        // this._childNodes is not updated, under the assumption that
        // ways are always downloaded with their child nodes.
    },

    // Updates calculated properties (parentWays, parentRels) for the specified change
    _updateCalculated: function(oldentity, entity, parentWays, parentRels) {

        parentWays = parentWays || this._parentWays;
        parentRels = parentRels || this._parentRels;

        var type = entity && entity.type || oldentity && oldentity.type,
            removed, added, ways, rels, i;


        if (type === 'way') {

            // Update parentWays
            if (oldentity && entity) {
                removed = _.difference(oldentity.nodes, entity.nodes);
                added = _.difference(entity.nodes, oldentity.nodes);
            } else if (oldentity) {
                removed = oldentity.nodes;
                added = [];
            } else if (entity) {
                removed = [];
                added = entity.nodes;
            }
            for (i = 0; i < removed.length; i++) {
                parentWays[removed[i]] = _.without(parentWays[removed[i]], oldentity.id);
            }
            for (i = 0; i < added.length; i++) {
                ways = _.without(parentWays[added[i]], entity.id);
                ways.push(entity.id);
                parentWays[added[i]] = ways;
            }

        } else if (type === 'relation') {

            // Update parentRels
            if (oldentity && entity) {
                removed = _.difference(oldentity.members, entity.members);
                added = _.difference(entity.members, oldentity);
            } else if (oldentity) {
                removed = oldentity.members;
                added = [];
            } else if (entity) {
                removed = [];
                added = entity.members;
            }
            for (i = 0; i < removed.length; i++) {
                parentRels[removed[i].id] = _.without(parentRels[removed[i].id], oldentity.id);
            }
            for (i = 0; i < added.length; i++) {
                rels = _.without(parentRels[added[i].id], entity.id);
                rels.push(entity.id);
                parentRels[added[i].id] = rels;
            }
        }
    },

    replace: function(entity) {
        if (this.entities[entity.id] === entity)
            return this;

        return this.update(function() {
            this._updateCalculated(this.entities[entity.id], entity);
            this.entities[entity.id] = entity;
        });
    },

    remove: function(entity) {
        return this.update(function() {
            this._updateCalculated(entity, undefined);
            this.entities[entity.id] = undefined;
        });
    },

    revert: function(id) {
        var baseEntity = this.base().entities[id],
            headEntity = this.entities[id];

        if (headEntity === baseEntity)
            return this;

        return this.update(function() {
            this._updateCalculated(headEntity, baseEntity);
            delete this.entities[id];
        });
    },

    update: function() {
        var graph = this.frozen ? iD.Graph(this, true) : this;

        for (var i = 0; i < arguments.length; i++) {
            arguments[i].call(graph, graph);
        }

        if (this.frozen) graph.frozen = true;

        return graph;
    },

    // Obliterates any existing entities
    load: function(entities) {
        var base = this.base();
        this.entities = Object.create(base.entities);

        for (var i in entities) {
            this.entities[i] = entities[i];
            this._updateCalculated(base.entities[i], this.entities[i]);
        }

        return this;
    }
};
iD.oneWayTags = {
    'aerialway': {
        'chair_lift': true,
        'mixed_lift': true,
        't-bar': true,
        'j-bar': true,
        'platter': true,
        'rope_tow': true,
        'magic_carpet': true,
        'yes': true
    },
    'highway': {
        'motorway': true,
        'motorway_link': true
    },
    'junction': {
        'roundabout': true
    },
    'man_made': {
        'piste:halfpipe': true
    },
    'piste:type': {
        'downhill': true,
        'sled': true,
        'yes': true
    },
    'waterway': {
        'river': true,
        'stream': true
    }
};
    module.exports =  {
        Difference: iD.Difference,
        Tree: iD.Tree,
        Entity: iD.Entity,
        Node: iD.Node,
        Way: iD.Way,
        Relation: iD.Relation,
        Graph: iD.Graph,
        geo: iD.geo,
        presets: iD.presets,
        locale: locale
};
