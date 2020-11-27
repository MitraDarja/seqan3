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
#include <seqan3/range/views/syncmer.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include <gtest/gtest.h>

#include "../iterator_test_template.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_shape;
using result_t = std::vector<size_t>;

inline static constexpr auto smer_view = seqan3::views::kmer_hash(seqan3::ungapped{3});
inline static constexpr auto gapped_smer_view = seqan3::views::kmer_hash(0b101_shape);
inline static constexpr auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{5});

//inline static constexpr auto syncmer_view = seqan3::views::syncmer(3,5);

using iterator_type = std::ranges::iterator_t< decltype(std::declval<seqan3::dna4_vector&>()
                                               | kmer_view
                                               | seqan3::views::syncmer( std::declval<seqan3::dna4_vector&>() | smer_view, 3, 0))>;

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true;

    seqan3::dna4_vector text{"AAGGCGT"_dna4};
    decltype(seqan3::views::kmer_hash(text, seqan3::ungapped{5})) vec = text | seqan3::views::kmer_hash(seqan3::ungapped{5});
    result_t expected_range{41, 166}; // AAGGC, AGGCG

    decltype(seqan3::views::syncmer(seqan3::views::kmer_hash(text, seqan3::ungapped{5}), seqan3::views::kmer_hash(text, seqan3::ungapped{2}), 2, 0)) test_range =
    seqan3::views::syncmer(vec, text | seqan3::views::kmer_hash(seqan3::ungapped{2}), 2, 0);
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
                                                std::list<seqan3::dna4> const,
                                                std::forward_list<seqan3::dna4>,
                                                std::forward_list<seqan3::dna4> const>;
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
    // ACGGC, CGGCG, ACGTT, CGTTT, GTTTA
    result_t result3{105, 422, 111, 447, 764};
    result_t result3_stop{105, 422};       // For stop at first T

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

    auto v = text | kmer_view | seqan3::views::syncmer(text | smer_view, 2, 0);
    compare_types(v);
}

TYPED_TEST(syncmer_view_properties_test, different_inputs_kmer_hash)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG
    // ACGTC, CGTCG, ACGTT, CGTTT, GTTTA
    result_t result{109, 438, 111, 447, 764};
    // TCGAC, GACGT
    result_t ungapped1{865, 539};
    // TCGAC, GACGT, TTTAG
    result_t gapped1{865, 539, 1010};
    EXPECT_RANGE_EQ(result, text | kmer_view | seqan3::views::syncmer(text | smer_view, 2, 0));
    EXPECT_RANGE_EQ(ungapped1, text | kmer_view | seqan3::views::syncmer(text | smer_view, 2, 1));
    EXPECT_RANGE_EQ(result, text | kmer_view | seqan3::views::syncmer(text | gapped_smer_view, 2, 0));
    EXPECT_RANGE_EQ(gapped1, text | kmer_view | seqan3::views::syncmer(text | gapped_smer_view, 2, 1));
}

TEST_F(syncmer_test, ungapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, text1 | kmer_view | seqan3::views::syncmer(text1 | smer_view, 2, 0) );
    EXPECT_RANGE_EQ(result1_short, text1_short | kmer_view | seqan3::views::syncmer(text1 | smer_view, 2, 0));
    auto empty_view = too_short_text | kmer_view | seqan3::views::syncmer(too_short_text | smer_view, 2, 0);
    EXPECT_TRUE(std::ranges::empty(empty_view));
    EXPECT_RANGE_EQ(result3, text3 | kmer_view | seqan3::views::syncmer(text3 | smer_view, 2, 0));
}

TEST_F(syncmer_test, gapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, text1 | kmer_view | seqan3::views::syncmer(text1 | gapped_smer_view, 2, 0));
    EXPECT_RANGE_EQ(result1_short, text1_short | kmer_view | seqan3::views::syncmer(text1 | gapped_smer_view, 2, 0));
    auto empty_view = too_short_text | kmer_view | seqan3::views::syncmer(too_short_text | gapped_smer_view, 2, 0);
    EXPECT_TRUE(std::ranges::empty(empty_view));
    EXPECT_RANGE_EQ(result3, text3 | kmer_view | seqan3::views::syncmer(text3 | gapped_smer_view, 2, 0));
}

TEST_F(syncmer_test, combinability)
{
    auto stop_at_t = seqan3::views::take_until([] (seqan3::dna4 const x) { return x == 'T'_dna4; });
    EXPECT_RANGE_EQ(result3_stop, text3 | stop_at_t | kmer_view | seqan3::views::syncmer(text3 | smer_view, 2, 0));
    EXPECT_RANGE_EQ(result3_stop, text3 | stop_at_t | kmer_view | seqan3::views::syncmer(text3 | gapped_smer_view, 2, 0));
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
    EXPECT_RANGE_EQ("AGGCGAGTTTAG"_dna4, text3 | seqan3::views::syncmer(text3_smers, 2, 0));
    //  CG
    // G
    // ACGCGAGTAG
}*/
