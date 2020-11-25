// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Mitra Darvish <mitra.darvish AT fu-berlin.de>
 * \brief Provides seqan3::views::syncmer_reverse_hash.
 */

#pragma once

#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/range/views/complement.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/syncmer.hpp>

#include <seqan3/range/views/zip.hpp>

#include <seqan3/core/debug_stream.hpp>

namespace seqan3
{
//!\brief strong_type for seed.
struct seed : seqan3::detail::strong_type<uint64_t, seed>
{
    using seqan3::detail::strong_type<uint64_t, seed>::strong_type;
};
} // namespace seqan3

namespace seqan3::detail
{
//!\brief seqan3::views::syncmer_reverse_hash's range adaptor object type (non-closure).
struct syncmer_reverse_hash_fn
{
    /*!\brief Store the shape and the window size and return a range adaptor closure object.
    * \param[in] shape       The seqan3::shape to use for hashing.
    * \param[in] window_size The windows size to use.
    * \throws std::invalid_argument if the size of the shape is greater than the `window_size`.
    * \returns               A range of converted elements.
    */
    constexpr auto operator()(shape const & kmer, shape const & window_size) const
    {
        return seqan3::detail::adaptor_from_functor{*this, kmer, window_size};
    }

    /*!\brief Store the shape, the window size and the seed and return a range adaptor closure object.
    * \param[in] shape       The seqan3::shape to use for hashing.
    * \param[in] window_size The size of the window.
    * \param[in] seed        The seed to use.
    * \throws std::invalid_argument if the size of the shape is greater than the `window_size`.
    * \returns               A range of converted elements.
    */
    constexpr auto operator()(shape const & kmer, shape const & window_size, seed const seed) const
    {
        return seqan3::detail::adaptor_from_functor{*this, kmer, window_size, seed};
    }

    /*!\brief Call the view's constructor with the underlying view, a seqan3::shape and a window size as argument.
     * \param[in] urange      The input range to process. Must model std::ranges::viewable_range and the reference type
     *                        of the range must model seqan3::semialphabet.
     * \param[in] shape       The seqan3::shape to use for hashing.
     * \param[in] window_size The size of the window.
     * \param[in] seed        The seed to use.
     * \throws std::invalid_argument if the size of the shape is greater than the `window_size`.
     * \returns               A range of converted elements.
     */
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange,
                              shape const & kmer,
                              shape const & window_size,
                              seed const seed = seed{0x8F3F73B5CF1C9ADE}) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
            "The range parameter to views::syncmer_reverse_hash cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng_t>,
            "The range parameter to views::syncmer_reverse_hash must model std::ranges::forward_range.");
        static_assert(semialphabet<std::ranges::range_reference_t<urng_t>>,
            "The range parameter to views::syncmer_reverse_hash must be over elements of seqan3::semialphabet.");

        if (kmer.size() < window_size.size())
            throw std::invalid_argument{"The size of the shape cannot be greater than the window size."};

        auto forward_strand = std::forward<urng_t>(urange) | seqan3::views::kmer_hash(kmer)
                                                           | std::views::transform([seed] (uint64_t i)
                                                                                  {return i ^ seed.get();});
        auto smers = std::forward<urng_t>(urange) | seqan3::views::kmer_hash(window_size)
                                                  | std::views::transform([seed] (uint64_t i)
                                                                                {return i ^ seed.get();});

        auto reverse_strand = std::forward<urng_t>(urange) | seqan3::views::complement
                                                           | std::views::reverse
                                                           | seqan3::views::kmer_hash(kmer)
                                                           | std::views::transform([seed] (uint64_t i)
                                                                                  {return i ^ seed.get();})
                                                           | std::views::reverse;

        auto reverse_smers = std::forward<urng_t>(urange) | seqan3::views::complement
                                                           | std::views::reverse
                                                           | seqan3::views::kmer_hash(window_size)
                                                           | std::views::transform([seed] (uint64_t i)
                                                                                  {return i ^ seed.get();})
                                                           | std::views::reverse;

        auto combined_strand = seqan3::views::zip(forward_strand, reverse_strand) | std::views::transform( [ ] (auto i) {return std::min(std::get<0>(i), std::get<1>(i));});

        return seqan3::detail::syncmer_view(combined_strand, smers, reverse_smers,  kmer.size() - window_size.size());

    }
};

} // namespace seqan3::detail

namespace seqan3::views
{

/*!\name Alphabet related views
 * \{
 */

/*!\brief                    Computes syncmers for a range with a given shape, window size and seed.
 * \tparam urng_t            The type of the range being processed. See below for requirements. [template parameter is
 *                           omitted in pipe notation]
 * \param[in] urange         The range being processed. [parameter is omitted in pipe notation]
 * \param[in] shape          The seqan3::shape that determines how to compute the hash value.
 * \param[in] window_size    The window size to use.
 * \param[in] seed           The seed used to skew the hash values. Default: 0x8F3F73B5CF1C9ADE.
 * \returns                  A range of `size_t` where each value is the syncmer of the resp. window.
 *                           See below for the properties of the returned range.
 * \ingroup views
 *
 * \details
 *
 * A sequence can be presented by a small number of k-mers (syncmers). For a given shape and window size all k-mers
 * are determined in the forward strand and the backward strand and only the lexicographically smallest k-mer is
 * returned for one window. This process is repeated over every possible window of a sequence. If consecutive windows
 * share a syncmer, it is saved only once.
 * For example, in the sequence "TAAAGTGCTAAA" for an ungapped shape of length 3 and a window size of 5 the first,
 * the second and the last window contain the same syncmer "AAA".
 * Because the syncmers of the first two consecutive windows also share the same position, storing this syncmer
 * twice is redundant and it is stored only once. The "AAA" syncmer of the last window on the other hand is stored,
 * since it is located at an other position than the previous "AAA" syncmer and hence storing the second
 * "AAA"-syncmer is not redundant but necessary.
 *
 * ### Non-lexicographical Minimisers by skewing the hash value with a seed
 *
 * It might happen that a syncmer changes only slightly when sliding the window over the sequence. For instance, when
 * a syncmer starts with a repetition of A’s, then in the next window it is highly likely that the syncmer will
 * start with a repetition of A’s as well. Because it is only one A shorter, depending on how long the repetition is
 * this might go on for multiple window shifts. Saving these only slightly different syncmer makes no sense because
 * they contain no new information about the underlying sequence.
 * Additionally, sequences with a repetition of A’s will be seen as more similar to each other than they actually are.
 * As [Marçais et al.](https://doi.org/10.1093/bioinformatics/btx235) have shown, randomizing the order of the k-mers
 * can solve this problem. Therefore, a random seed is used to XOR all k-mers, thereby randomizing the
 * order. The user can change the seed to any other value he or she thinks is useful. A seed of 0 is returning the
 * lexicographical order.
 *
 * \sa seqan3::views::syncmer_view
 *
 * \attention
 * Be aware of the requirements of the seqan3::views::kmer_hash view.
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)   | `rrng_t` (returned range type)   |
 * |----------------------------------|:----------------------------------:|:--------------------------------:|
 * | std::ranges::input_range         | *required*                         | *preserved*                      |
 * | std::ranges::forward_range       | *required*                         | *preserved*                      |
 * | std::ranges::bidirectional_range |                                    | *lost*                           |
 * | std::ranges::random_access_range |                                    | *lost*                           |
 * | std::ranges::contiguous_range    |                                    | *lost*                           |
 * |                                  |                                    |                                  |
 * | std::ranges::viewable_range      | *required*                         | *guaranteed*                     |
 * | std::ranges::view                |                                    | *guaranteed*                     |
 * | std::ranges::sized_range         |                                    | *lost*                           |
 * | std::ranges::common_range        |                                    | *lost*                           |
 * | std::ranges::output_range        |                                    | *lost*                           |
 * | seqan3::const_iterable_range     |                                    | *preserved*                      |
 * |                                  |                                    |                                  |
 * | std::ranges::range_reference_t   | seqan3::semialphabet               | std::size_t                      |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * \include test/snippet/range/views/syncmer_reverse_hash.cpp
 *
 * \hideinitializer
 */
inline constexpr auto syncmer_reverse_hash = detail::syncmer_reverse_hash_fn{};

//!\}

} // namespace seqan3::views
