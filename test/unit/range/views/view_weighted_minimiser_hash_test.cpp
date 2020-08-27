// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <forward_list>
#include <list>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>
#include <seqan3/range/views/weighted_minimiser_hash.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include <gtest/gtest.h>

#include "../iterator_test_template.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_shape;
using result_t = std::vector<size_t>;



static constexpr seqan3::shape ungapped_shape = seqan3::ungapped{4};
static constexpr seqan3::shape gapped_shape = 0b1001_shape;

seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> bloomfilter {seqan3::bin_count{1u},
                                                                                seqan3::bin_size{100u},
                                                                                seqan3::hash_function_count{1}};

using iterator_type = std::ranges::iterator_t<decltype(seqan3::views::weighted_minimiser_hash(
                                                       std::declval<seqan3::dna4_vector&>(),
                                                       ungapped_shape,
                                                       seqan3::window_size{8},
                                                       bloomfilter,
                                                       seqan3::seed{0}))>;

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = false;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    result_t expected_range{26, 97, 27, 6};

    using test_range_t = decltype(seqan3::views::weighted_minimiser_hash(text,
                                                                         ungapped_shape,
                                                                         seqan3::window_size{8},
                                                                         bloomfilter,
                                                                         seqan3::seed{0}));
    test_range_t test_range = seqan3::views::weighted_minimiser_hash(text,
                                                                     ungapped_shape,
                                                                     seqan3::window_size{8},
                                                                     bloomfilter,
                                                                     seqan3::seed{0});
};

using test_type = ::testing::Types<iterator_type>;
INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_type, );

template <typename T>
class weighted_minimiser_hash_properties_test: public ::testing::Test { };

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                seqan3::bitcompressed_vector<seqan3::dna4>,
                                                seqan3::bitcompressed_vector<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const>;

TYPED_TEST_SUITE(weighted_minimiser_hash_properties_test, underlying_range_types, );
class weighted_minimiser_hash_test : public ::testing::Test
{
protected:
    std::vector<seqan3::dna4> text1{"AAAAAAAAAAAAAAAAAAA"_dna4};
    std::vector<seqan3::dna4> text1_short{"AAAAAA"_dna4};
    result_t result1{0, 0, 0}; // Same for ungapped and gapped
    result_t ungapped_default_seed{0x8F3F73B5CF1C9A21, 0x8F3F73B5CF1C9A21, 0x8F3F73B5CF1C9A21};
    result_t gapped_default_seed{0x8F3F73B5CF1C9AD1, 0x8F3F73B5CF1C9AD1, 0x8F3F73B5CF1C9AD1};

    std::vector<seqan3::dna4> text2{"AC"_dna4};
    result_t result2{};

    std::vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4};
    result_t ungapped3{26, 101, 27, 6};    // ACGG, CGCC, ACGT, aacg, aaac
    result_t ungapped_stop_at_t3{26, 101}; // ACGG, CGCC
    result_t gapped3{2, 5, 3, 2};          // A--G, C--C, A--T, A--G "-" for gap
    result_t gapped_stop_at_t3{2, 5};      // A--G, C--C "-" for gap
};


TYPED_TEST(weighted_minimiser_hash_properties_test, different_input_ranges)
{

    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                   'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG
    result_t ungapped{27, 109, 27, 6}; // ACGT, CGTC, ACGT, aacg
    result_t gapped{3, 5, 3, 2};      // A--T, C--C, A--T, a--g - "-" for gap

    bloomfilter.emplace(97, seqan3::bin_index{0u}); // CGAC
    bloomfilter.emplace(1, seqan3::bin_index{0u});  // aaac

    EXPECT_RANGE_EQ(ungapped, seqan3::views::weighted_minimiser_hash(text,
                                                                     ungapped_shape,
                                                                     seqan3::window_size{8},
                                                                     bloomfilter,
                                                                     seqan3::seed{0}));
    EXPECT_RANGE_EQ(gapped, seqan3::views::weighted_minimiser_hash(text,
                                                                   gapped_shape,
                                                                   seqan3::window_size{8},
                                                                   bloomfilter,
                                                                   seqan3::seed{0}));
}

TEST_F(weighted_minimiser_hash_test, ungapped)
{
    bloomfilter.emplace(97, seqan3::bin_index{0u}); // CGAC
    bloomfilter.emplace(1, seqan3::bin_index{0u});  // aaac
    EXPECT_RANGE_EQ(result1, seqan3::views::weighted_minimiser_hash(text1,
                                                                    ungapped_shape,
                                                                    seqan3::window_size{8},
                                                                    bloomfilter,
                                                                    seqan3::seed{0}));
    EXPECT_RANGE_EQ(result2, seqan3::views::weighted_minimiser_hash(text2,
                                                                    ungapped_shape,
                                                                    seqan3::window_size{8},
                                                                    bloomfilter,
                                                                    seqan3::seed{0}));
    EXPECT_RANGE_EQ(ungapped3, seqan3::views::weighted_minimiser_hash(text3,
                                                                      ungapped_shape,
                                                                      seqan3::window_size{8},
                                                                      bloomfilter,
                                                                      seqan3::seed{0}));

    auto stop_at_t = seqan3::views::take_until([] (seqan3::dna4 const x) { return x == 'T'_dna4; });
    EXPECT_RANGE_EQ(ungapped_stop_at_t3, seqan3::views::weighted_minimiser_hash(text3 | stop_at_t,
                                                                                ungapped_shape,
                                                                                seqan3::window_size{8},
                                                                                bloomfilter,
                                                                                seqan3::seed{0}));
}

TEST_F(weighted_minimiser_hash_test, gapped)
{
    bloomfilter.emplace(97, seqan3::bin_index{0u}); // CGAC
    bloomfilter.emplace(1, seqan3::bin_index{0u});  // aaac
    EXPECT_RANGE_EQ(result1, seqan3::views::weighted_minimiser_hash(text1,
                                                                    gapped_shape,
                                                                    seqan3::window_size{8},
                                                                    bloomfilter,
                                                                    seqan3::seed{0}));
    EXPECT_RANGE_EQ(result2, seqan3::views::weighted_minimiser_hash(text2,
                                                                    gapped_shape,
                                                                    seqan3::window_size{8},
                                                                    bloomfilter,
                                                                    seqan3::seed{0}));
    EXPECT_RANGE_EQ(gapped3, seqan3::views::weighted_minimiser_hash(text3,
                                                                      gapped_shape,
                                                                      seqan3::window_size{8},
                                                                      bloomfilter,
                                                                      seqan3::seed{0}));

    auto stop_at_t = seqan3::views::take_until([] (seqan3::dna4 const x) { return x == 'T'_dna4; });
    EXPECT_RANGE_EQ(gapped_stop_at_t3, seqan3::views::weighted_minimiser_hash(text3 | stop_at_t,
                                                                                gapped_shape,
                                                                                seqan3::window_size{8},
                                                                                bloomfilter,
                                                                                seqan3::seed{0}));
}


TEST_F(weighted_minimiser_hash_test, seed)
{
    EXPECT_RANGE_EQ(ungapped_default_seed,
                    seqan3::views::weighted_minimiser_hash(text1, ungapped_shape, seqan3::window_size{8}, bloomfilter));
    EXPECT_RANGE_EQ(gapped_default_seed, seqan3::views::weighted_minimiser_hash(text1,
                                                                                gapped_shape,
                                                                                seqan3::window_size{8},
                                                                                bloomfilter));
}

TEST_F(weighted_minimiser_hash_test, shape_bigger_than_window)
{
    EXPECT_THROW(seqan3::views::weighted_minimiser_hash(text1, ungapped_shape, seqan3::window_size{3}, bloomfilter),
                                                        std::invalid_argument);
    EXPECT_THROW(seqan3::views::weighted_minimiser_hash(text1, gapped_shape, seqan3::window_size{3}, bloomfilter),
                                                        std::invalid_argument);
}