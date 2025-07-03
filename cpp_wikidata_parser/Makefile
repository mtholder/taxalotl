includedir=picojson

parse_full_wikidata_dump: picojson/picojson.h parse_full_wikidata_dump.cpp
	$(CXX) parse_full_wikidata_dump.cpp -DPICOJSON_USE_INT64 -O4 -Wall -o parse_full_wikidata_dump


clean:
	rm -f parse_full_wikidata_dump


