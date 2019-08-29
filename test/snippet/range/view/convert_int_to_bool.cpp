#include <vector>

#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/std/ranges>

int main()
{
    // convert from int to bool
    std::vector<int>  vec{7, 5, 0, 5, 0, 0, 4, 8, -3};

    // pipe notation
    auto v = vec | seqan3::view::convert<bool>; // == [1, 1, 0, 1, 0, 0, 1, 1, 1];

    // function notation and immediate conversion to vector again
    auto v2 = seqan3::view::convert<bool>(vec) | std::ranges::to<std::vector<bool>>;

    // combinability
    auto v3 = vec | seqan3::view::convert<bool> | std::view::reverse; // == [1, 1, 1, 0, 0, 1, 0, 1, 1];
}