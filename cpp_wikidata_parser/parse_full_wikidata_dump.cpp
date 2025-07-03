// Modified by MTH from the examples/streaming.cc from
//    the picojson project at https://github.com/kazuho/picojson.git
/*
 * Copyright 2009-2010 Cybozu Labs, Inc.
 * Copyright 2011-2014 Kazuho Oku
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#include <iostream>
#include <iterator>
#include <string>
#include "picojson/picojson.h"

using std::endl;
using std::cerr;
using std::cout;

class root_context {
  public:
    bool parse_array_start() {
        cerr << "parse_array_start" << endl;
        return true; // only allow array as root
    }

    template <typename Iter> bool parse_array_item(picojson::input<Iter> &in, 
                                                   size_t) {
        cerr << "parse_array_item" << endl;
        ++el_count;
        picojson::value item;
        // parse the array item
        picojson::default_parse_context ctx(&item);
        if (!picojson::_parse(ctx, in)) {
          return false;
        }
        // assert that the array item is a hash
        if (!item.is<picojson::object>()) {
          return false;
        }
        std::string idstr = item.get("id").to_str();
        
        const picojson::value & labelsval = item.get("labels");
        const picojson::value & enlabel = labelsval.get("en");
        const picojson::value & english_name_val = enlabel.get("value");
        std::string english_name_str = english_name_val.to_str();
        cout << "id=" << idstr << "\t" << english_name_str <<  endl;
        return true;
    }

    bool parse_object_stop() {
        cerr << "parse_object_stop" << endl;
        return true;
    }

    bool parse_array_stop(size_t) {
        cerr << "parse_array_stop" << endl;
        return true;
    }

    size_t get_el_count() const {
        return this->el_count;
    }

    bool set_null() {
        cerr << "set_null" << endl;
        return true;
    }

    bool set_bool(bool v) {
        cerr << "set_bool(" << v << ")" << endl;
        return true;
    }
#   ifdef PICOJSON_USE_INT64
        bool set_int64(int64_t v) {
            cerr << "set_int64(" << v << ")" << endl;
            return true;
        }
#   endif
    bool set_number(double v) {
        cerr << "set_number(" << v << ")" << endl;
        return true;
    }

    template <typename Iter> bool parse_string(picojson::input<Iter> &) {
        cerr << "parse_string" << endl;
        return true;
    }
  
    bool parse_object_start() {
        cerr << "parse_object_start" << endl;
        return true;
    }

    template <typename Iter> bool parse_object_item(picojson::input<Iter> &, 
    const std::string &) {
        cerr << "parse_object_item" << endl;
        return true;
    }
  protected:
    size_t el_count = 0;
};


int main(void) {
  root_context ctx;
  std::string err;

  picojson::_parse(ctx,
                   std::istream_iterator<char>(std::cin),
                   std::istream_iterator<char>(), &err);

  if (!err.empty()) {
    std::cerr  << "\nError above.\n" << std::endl;
    return 1;
  }
  std::cout << ctx.get_el_count() << " elements read." << std::endl;

  return 0;
}
