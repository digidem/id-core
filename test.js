var chai = require('chai')
var expect = global.expect = chai.expect;
var _ = global._ = require('lodash');

chai.use(function (chai, utils) {
    var flag = utils.flag;

    chai.Assertion.addMethod('classed', function (className) {
        this.assert(
            flag(this, 'object').classed(className)
            , 'expected #{this} to be classed #{exp}'
            , 'expected #{this} not to be classed #{exp}'
            , className
        );
    });
});

// Tests depend on the iD editor default OSM presets
var presets = {
    presets: require('iD/data/presets/presets'),
    defaults: require('iD/data/presets/defaults'),
    categories: require('iD/data/presets/categories'),
    fields: require('iD/data/presets/fields')
}

var iD = global.iD = require('./')(presets);

iD.data = {
    presets: presets
};

// #search tests on iD.presets.Collection depend on English translations
iD.locale.en = require('iD/dist/locales/en.json');
iD.locale.current('en');

require('iD/test/spec/core/difference');
require('iD/test/spec/core/entity');
require('iD/test/spec/core/graph');
require('iD/test/spec/core/node');
require('iD/test/spec/core/relation');
require('iD/test/spec/core/tree');
require('iD/test/spec/core/way');

require('iD/test/spec/presets');
require('iD/test/spec/presets/category');
require('iD/test/spec/presets/collection');
require('iD/test/spec/presets/preset');

require('iD/test/spec/geo');
require('iD/test/spec/geo/extent');
require('iD/test/spec/geo/intersection');
require('iD/test/spec/geo/multipolygon');
