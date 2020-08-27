// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Mitra Darvish <mitra.darvish AT fu-berlin.de>
 * \brief Provides seqan3::views::weighted_minimiser_hash.
 */

#pragma once

#include <seqan3/range/views/minimiser_hash.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/range/views/zip.hpp>

namespace seqan3::detail
{
//!\brief seqan3::views::weighted_minimiser_hash's range adaptor object type (non-closure).
struct weighted_minimiser_hash_fn
{
    /*!\brief Store the shape and the window size and return a range adaptor closure object.
    * \param[in] shape       The seqan3::shape to use for hashing.
    * \param[in] window_size The windows size to use.
    * \param[in] bloomfilter The bloomfilter with hash values that should be down weighted.
    * \throws std::invalid_argument if the size of the shape is greater than the `window_size`.
    * \returns               A range of converted elements.
    */
    /*template <class IBFType>
    constexpr auto operator()(shape const & shape, window_size const window_size,  IBFType const & bloomfilter) const
    {
        return seqan3::detail::adaptor_from_functor{*this, shape, window_size, bloomfilter};
    }*/

    /*!\brief Store the shape, the window size and the seed and return a range adaptor closure object.
    * \param[in] shape       The seqan3::shape to use for hashing.
    * \param[in] window_size The size of the window.
    * \param[in] bloomfilter The bloomfilter with hash values that should be down weighted.
    * \param[in] seed        The seed to use.
    * \throws std::invalid_argument if the size of the shape is greater than the `window_size`.
    * \returns               A range of converted elements.
    */
    /*template <class IBFType>
    constexpr auto operator()(shape const & shape,
                              window_size const window_size,
                              IBFType const & bloomfilter,
                              seed const seed) const
    {
        return seqan3::detail::adaptor_from_functor{*this, shape, window_size, bloomfilter, seed};
    }*/

    /*!\brief Call the view's constructor with the underlying view, a seqan3::shape and a window size as argument.
     * \param[in] urange      The input range to process. Must model std::ranges::viewable_range and the reference type
     *                        of the range must model seqan3::semialphabet.
     * \param[in] shape       The seqan3::shape to use for hashing.
     * \param[in] window_size The size of the window.
     * \param[in] bloomfilter The bloomfilter with hash values that should be down weighted.
     * \param[in] seed        The seed to use.
     * \throws std::invalid_argument if the size of the shape is greater than the `window_size`.
     * \returns               A range of converted elements.
     */
    template <std::ranges::range urng_t, class IBFType>
    constexpr auto operator()(urng_t && urange,
                              shape const & shape,
                              window_size const window_size,
                              IBFType const & bloomfilter,
                              seed const seed = seed{0x8F3F73B5CF1C9ADE}) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
            "The range parameter to views::weighted_weighted_minimiser_hash cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng_t>,
            "The range parameter to views::weighted_minimiser_hash must model std::ranges::forward_range.");
        static_assert(semialphabet<std::ranges::range_reference_t<urng_t>>,
            "The range parameter to views::weighted_minimiser_hash must be over elements of seqan3::semialphabet.");

        if (shape.size() > window_size.get())
            throw std::invalid_argument{"The size of the shape cannot be greater than the window size."};

        auto forward_strand = std::forward<urng_t>(urange) | seqan3::views::kmer_hash(shape)
                                                           | std::views::transform([seed] (uint64_t i)
                                                                                  {return i ^ seed.get();});

        auto reverse_strand = std::forward<urng_t>(urange) | seqan3::views::complement
                                                           | std::views::reverse
                                                           | seqan3::views::kmer_hash(shape)
                                                           | std::views::transform([seed] (uint64_t i)
                                                                                  {return i ^ seed.get();})
                                                           | std::views::reverse;



        auto both = seqan3::views::zip(forward_strand, reverse_strand) | std::views::transform( [bloomfilter] (auto i)
                    {
                        auto agent = bloomfilter.membership_agent();
                        if ((agent.bulk_contains(std::get<0>(i))[0] > 0) | (agent.bulk_contains(std::get<1>(i))[0] > 0))
                            return std::max(std::get<0>(i), std::get<1>(i));
                        else
                            return std::min(std::get<0>(i), std::get<1>(i));
                    });

        return seqan3::detail::minimiser_view(both, window_size.get() - shape.size() + 1);
    }
};

} // namespace seqan3::detail

namespace seqan3::views
{

/*!\name Alphabet related views
 * \{
 */

/*!\brief                    Computes minimisers for a range with a given shape, window size, seed and a known list of
 *                           unfavourable k-mers.
 * \tparam urng_t            The type of the range being processed. See below for requirements. [template parameter is
 *                           omitted in pipe notation]
 * \param[in] urange         The range being processed. [parameter is omitted in pipe notation]
 * \param[in] shape          The seqan3::shape that determines how to compute the hash value.
 * \param[in] window_size    The window size to use.
 * \param[in] bloomfilter    The bloomfilter with hash values that should be down weighted.
 * \param[in] seed           The seed used to skew the hash values. Default: 0x8F3F73B5CF1C9ADE.
 * \returns                  A range of `size_t` where each value is the minimiser of the resp. window.
 *                           See below for the properties of the returned range.
 * \ingroup views
 *
 * \details
 *
 * This view is really similar to seqan3::views::minimiser_hash, see there for an explaination of minimisers and an
 * explanation for the usefulness of seeds. This view adds another parameter, a bloomfilter storing k-mers which are
 * less favourable and therefore should be less likely to be picked as minimisers. The general idea was presented by
 * [Jain et al.](https://www.biorxiv.org/content/10.1101/2020.02.11.943241v1.full.pdf),
 * the implementation here uses another hash function though and therefore another weighting system.
 * If a k-mer or its reverse complement is found in the list of unfavourable k-mers, then the maximum hash value
 * between the hash value of the k-mer and the hash value of the reverse complement is considered instead of its
 * minimum. Thereby making it less likely that the value is considered a minimiser. Note: Due to this definition there
 * is no possiblity to make k-mers less likely when their reverse complement is identical.
 *
 * \sa seqan3::views::minimiser_view
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
 *
 * \hideinitializer
 */
inline constexpr auto weighted_minimiser_hash = detail::weighted_minimiser_hash_fn{};

//!\}

} // namespace seqan3::views
