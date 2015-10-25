# See the README for installation instructions.

all: \
	lib/d3.js \
	index.js

clean:
	rm -f lib/d3.js index.js

D3_FILES = \
	node_modules/d3/src/start.js \
	node_modules/d3/src/arrays/values.js \
	node_modules/d3/src/geo/area.js \
	node_modules/d3/src/end.js

lib/d3.js: $(D3_FILES)
	node_modules/.bin/smash $(D3_FILES) > $@

ID_CORE_FILES = \
	start.js \
  node_modules/iD/js/id/geo.js \
  node_modules/iD/js/id/geo/extent.js \
  node_modules/iD/js/id/geo/intersection.js \
  node_modules/iD/js/id/geo/multipolygon.js \
  node_modules/iD/js/id/geo/raw_mercator.js \
  node_modules/iD/js/id/actions/reverse.js \
  node_modules/iD/js/id/util.js \
  node_modules/iD/js/id/presets.js \
  node_modules/iD/js/id/presets/category.js \
  node_modules/iD/js/id/presets/collection.js \
  node_modules/iD/js/id/presets/field.js \
  node_modules/iD/js/id/presets/preset.js \
  node_modules/iD/js/id/core/difference.js \
  node_modules/iD/js/id/core/tree.js \
  node_modules/iD/js/id/core/entity.js \
  node_modules/iD/js/id/core/node.js \
  node_modules/iD/js/id/core/way.js \
  node_modules/iD/js/id/core/relation.js \
  node_modules/iD/js/id/core/graph.js \
  node_modules/iD/js/id/core/oneway_tags.js \
  end.js

index.js: $(ID_CORE_FILES)
	node_modules/.bin/smash $(ID_CORE_FILES) > $@
