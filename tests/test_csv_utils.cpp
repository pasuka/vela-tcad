#include "vela/io/CsvUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <stdexcept>

TEST_CASE("CSV utility splits and trims unquoted rows", "[csv]")
{
    const auto columns = vela::splitCsvLine(" node_id, donors_cm3, acceptors_cm3\r");

    REQUIRE(columns.size() == 3);
    REQUIRE(columns[0] == "node_id");
    REQUIRE(columns[1] == "donors_cm3");
    REQUIRE(columns[2] == "acceptors_cm3");
}

TEST_CASE("CSV utility rejects quoted fields", "[csv]")
{
    REQUIRE_THROWS_AS(vela::splitCsvLine("node_id,\"donors_cm3\",acceptors_cm3"),
                      std::runtime_error);
}
