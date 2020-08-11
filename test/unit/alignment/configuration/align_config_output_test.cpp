// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/configuration/align_config_output.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/type_traits/basic.hpp>

#include "../../core/algorithm/pipeable_config_element_test_template.hpp"

template <typename test_t>
struct align_cfg_output_test : public ::testing::Test
{};

using test_types = ::testing::Types<seqan3::align_cfg::output_score_tag,
                                    seqan3::align_cfg::output_end_position_tag,
                                    seqan3::align_cfg::output_begin_position_tag>;

INSTANTIATE_TYPED_TEST_SUITE_P(output, pipeable_config_element_test, test_types, );

TEST(align_config_output, score)
{
    EXPECT_TRUE((std::same_as<seqan3::remove_cvref_t<decltype(seqan3::align_cfg::output_score)>,
                              seqan3::align_cfg::output_score_tag>));
}

TEST(align_config_output, end_position)
{
    EXPECT_TRUE((std::same_as<seqan3::remove_cvref_t<decltype(seqan3::align_cfg::output_end_position)>,
                              seqan3::align_cfg::output_end_position_tag>));
}

TEST(align_config_output, begin_position)
{
    EXPECT_TRUE((std::same_as<seqan3::remove_cvref_t<decltype(seqan3::align_cfg::output_begin_position)>,
                              seqan3::align_cfg::output_begin_position_tag>));
}

TEST(align_config_output, combine_outputs)
{
    seqan3::configuration cfg = seqan3::align_cfg::output_score |
                                seqan3::align_cfg::output_end_position |
                                seqan3::align_cfg::output_begin_position;

    EXPECT_TRUE(cfg.exists<seqan3::align_cfg::output_score_tag>());
    EXPECT_TRUE(cfg.exists<seqan3::align_cfg::output_end_position_tag>());
    EXPECT_TRUE(cfg.exists<seqan3::align_cfg::output_begin_position_tag>());
}
