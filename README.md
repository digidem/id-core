# iD-core - core data types from [iDEditor](https://github.com/openstreetmap/iD)

[![build status](https://img.shields.io/travis/digidem/id-core/master.svg)](https://travis-ci.org/digidem/id-core)
[![npm version](https://img.shields.io/npm/v/id-core.svg)](https://www.npmjs.com/package/id-core)

**ALPHA:** see [TODO](#TODO)

## Overview

This packages up [iD Editor](https://github.com/openstreetmap/iD) [core OSM data types and Graph](https://github.com/openstreetmap/iD/blob/master/ARCHITECTURE.md#core) as an npm package to `require()` in your own projects. Included:

* [iD.Difference](https://github.com/openstreetmap/iD/blob/master/js/id/core/difference.js)
* [iD.Tree](https://github.com/openstreetmap/iD/blob/master/js/id/core/tree.js)
* [iD.Entity](https://github.com/openstreetmap/iD/blob/master/js/id/core/entity.js)
* [iD.Node](https://github.com/openstreetmap/iD/blob/master/js/id/core/node.js)
* [iD.Way](https://github.com/openstreetmap/iD/blob/master/js/id/core/way.js)
* [iD.Relation](https://github.com/openstreetmap/iD/blob/master/js/id/core/relation.js)
* [iD.Graph](https://github.com/openstreetmap/iD/blob/master/js/id/core/graph.js)
* [iD.geo](https://github.com/openstreetmap/iD/blob/master/js/id/geo.js)
* [iD.presets](https://github.com/openstreetmap/iD/blob/master/js/id/presets.js)
* [iD.locale](https://github.com/openstreetmap/iD/blob/master/js/lib/locale.js)

## Usage

Needs to be instantiated with [presets](https://github.com/openstreetmap/iD/tree/master/data/presets) **subject to change see [TODO](#TODO).


```sh
npm install id-core
```

```js
var presets = {
    categories: {},
    defaults: {},
    fields: {},
    presets: {}
};

var iDCore = require('id-core');

var myNode = new iDCore.Node();
```

## Why

According to https://github.com/openstreetmap/iD/blob/master/ARCHITECTURE.md#core

> [iD] eventually aims to be a reusable, modular library to kickstart other
> JavaScript-based tools for OpenStreetMap.

The OSM data model is complex and hard to implement. iD is not published on npm and importing the whole iD project is excessive for a JavaScript based tool for OpenStreetMap. At [Digital Democracy](http://www.digital-democracy.org/) we are building tools on top of OSM, and borrowing from iD gives us a head start.

## How

iD does not use a commonJS module structure, so it's not as simple as `require`ing what is needed. We use [Smash](https://github.com/mbostock/smash) to concatenate just what is needed from d3 and iD editor to make things work. To rebuild from iD source files:

```sh
make clean && make
```

## Tests

```sh
npm install
npm test
```

Uses tests directly from iD to test exported objects.

## TODO

* `iD-core` needs to be instantiated with presets, since the core types depend on them being set. We currently export the uninstantiated object though, which is needed for tests to work, but that means we have no access to the instantiated presets object.
* Needs documentation beyond https://github.com/openstreetmap/iD/blob/master/ARCHITECTURE.md#core
* `iD.locale` has no tests yet. It is adapted from [iD/js/lib/locale.js](https://github.com/openstreetmap/iD/blob/master/js/lib/locale.js) which polutes the global namespace, which we don't want. We can writed our own or wait for [iD/#890](https://github.com/openstreetmap/iD/issues/890).

## License

iD-core is available under the [WTFPL](http://sam.zoy.org/wtfpl/), though obviously,
if you want to dual-license any contributions that's cool.
It includes [d3js](http://d3js.org/), which BSD-licensed.
