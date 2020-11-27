// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <forward_list>
#include <list>
#include <type_traits>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/detail/debug_stream_alphabet.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>
#include <seqan3/range/views/complement.hpp>
#include <seqan3/range/views/drop.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/syncmer_reverse_hash.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include <gtest/gtest.h>

#include "../iterator_test_template.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_shape;
using result_t = std::vector<size_t>;


using iterator_type = std::ranges::iterator_t<decltype(std::declval<seqan3::dna4_vector&>()
                                                      | seqan3::views::syncmer_reverse_hash(seqan3::ungapped{5},
                                                                                    seqan3::ungapped{2},
                                                                                    0,
                                                                                    seqan3::seed{0}))>;

static constexpr seqan3::shape kmers = seqan3::ungapped{5};
static constexpr seqan3::shape ungapped_shape = seqan3::ungapped{3};
static constexpr seqan3::shape gapped_shape = 0b101_shape;
static constexpr auto ungapped_view = seqan3::views::syncmer_reverse_hash(kmers,
                                                                  ungapped_shape,
                                                                  0,
                                                                  seqan3::seed{0});
static constexpr auto gapped_view = seqan3::views::syncmer_reverse_hash(kmers,
                                                                gapped_shape,
                                                                0,
                                                                seqan3::seed{0});

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true; //TODO: Check, why it does not work with true!

    seqan3::dna4_vector text{"GGCAAGT"_dna4};
    result_t expected_range{505, 126}; // cttgc, acttg

    using test_range_t = decltype(text | seqan3::views::syncmer_reverse_hash(kmers,
                                                                      seqan3::ungapped{2},
                                                                      0,
                                                                      seqan3::seed{0}));
    test_range_t test_range = text |  seqan3::views::syncmer_reverse_hash(kmers,
                                                                      seqan3::ungapped{2},
                                                                      0,
                                                                      seqan3::seed{0});
};

using test_types = ::testing::Types<iterator_type>;
INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_types, );

template <typename T>
class syncmer_view_properties_test: public ::testing::Test { };

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                seqan3::bitcompressed_vector<seqan3::dna4>,
                                                seqan3::bitcompressed_vector<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const>;
                                                //std::forward_list<seqan3::dna4>, // TODO: Add them, at the moment segfault
                                                //std::forward_list<seqan3::dna4> const>;
TYPED_TEST_SUITE(syncmer_view_properties_test, underlying_range_types, );

class syncmer_test : public ::testing::Test
{
protected:
    std::vector<seqan3::dna4> text1{"AAAAAAAAAAAAAAAAAAA"_dna4};
    std::vector<seqan3::dna4> text1_short{"AAAAAA"_dna4};
    result_t result1{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    result_t result1_short{0, 0};

    std::vector<seqan3::dna4> too_short_text{"AC"_dna4};

    std::vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4};
    // ACGGC, CGGCG, cgtcg, CGACG, acgtc, aacgt, aaacg, GTTTA
    result_t result3{105, 406, 390, 109, 27, 6, 764};
    result_t result3_stop{105, 406, 390};       // For stop at first T

};

template <typename adaptor_t>
void compare_types(adaptor_t v)
{
    EXPECT_TRUE(std::ranges::input_range<decltype(v)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v)>);
    EXPECT_TRUE(std::ranges::view<decltype(v)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v), size_t>));
}

TYPED_TEST(syncmer_view_properties_test, concepts)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                   'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG

    auto v = text | ungapped_view;
    compare_types(v);
    auto v2 = text | gapped_view;
    compare_types(v2);

}


TYPED_TEST(syncmer_view_properties_test, different_inputs_kmer_hash)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG
    // ACGTC, CGTCG, CGTCG, acgtc, aacgt, aaacg, GTTTA
    result_t result{109, 390, 390, 109, 27, 6, 764};
    EXPECT_RANGE_EQ(result, text | ungapped_view);
    EXPECT_RANGE_EQ(result, text | gapped_view);
}

TEST_F(syncmer_test, ungapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, text1 | ungapped_view );
    EXPECT_RANGE_EQ(result1_short, text1_short | ungapped_view);
    auto empty_view = too_short_text | ungapped_view;
    EXPECT_TRUE(std::ranges::empty(empty_view));
    EXPECT_RANGE_EQ(result3, text3 | ungapped_view);
}

TEST_F(syncmer_test, gapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, text1 | gapped_view);
    EXPECT_RANGE_EQ(result1_short, text1_short | gapped_view);
    auto empty_view = too_short_text | gapped_view;
    EXPECT_TRUE(std::ranges::empty(empty_view));
    EXPECT_RANGE_EQ(result3, text3 | gapped_view);
}

TEST_F(syncmer_test, combinability)
{
    auto stop_at_t = seqan3::views::take_until([] (seqan3::dna4 const x) { return x == 'T'_dna4; });
    EXPECT_RANGE_EQ(result3_stop, text3 | stop_at_t | ungapped_view);
    EXPECT_RANGE_EQ(result3_stop, text3 | stop_at_t | gapped_view);
}
//TODO: This leads to a segfault, why?
/*
TEST_F(syncmer_test, non_arithmetic_value)
{
    // just compute the syncmers directly on the alphabet
    // ACGGCGACGTTTAG
    std::vector<seqan3::dna4> text3_smers{"AAAGCGGCGCGCGGGGCGCGGCTTAT"_dna4};
    //seqan3::debug_stream << text3_smers;
    seqan3::debug_stream << *std::ranges::begin(text3_smers) << "\n";
    seqan3::debug_stream <<  text3_smers << "\n";
    EXPECT_RANGE_EQ("AGGCGAGTTTAG"_dna4, text3 | seqan3::views::syncmer(text3_smers, 2));
    //  CG
    // G
    // ACGCGAGTAG
}*/
