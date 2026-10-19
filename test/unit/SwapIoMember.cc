/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "ijedi/Io/IoBase.h"

namespace {

CASE("swapIoMember/replaces_pattern_with_padded_member") {
  // Config as OOPS passes it for member 2: member sits beside the io block
  eckit::LocalConfiguration io;
  io.set("filetype", "unstructured");
  io.set("filepath", "./inc/mem%{member}%/increment");
  eckit::LocalConfiguration conf;
  conf.set("io", io);
  conf.set("member", 2);

  const eckit::LocalConfiguration out(ijedi::swapIoMember(conf), "io");
  EXPECT(out.getString("filepath") == "./inc/mem002/increment");
  EXPECT(out.getString("filetype") == "unstructured");
}

}  // namespace

int main(int argc, char ** argv) { return eckit::testing::run_tests(argc, argv); }
