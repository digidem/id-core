module.exports = function(_) {
    var presets = iD.presets()
    presets.load(_);
    iD.areaKeys = presets.areaKeys()

    return {
        Difference: iD.Difference,
        Tree: iD.Tree,
        Entity: iD.Entity,
        Node: iD.Node,
        Way: iD.Way,
        Relation: iD.Relation,
        Graph: iD.Graph,
        geo: iD.geo,
        presets: iD.presets,
        areaKeys: iD.areaKeys,
        locale: locale
    };
};
