// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::take.
 */

#pragma once

#include <range/v3/algorithm/copy.hpp>

#include <seqan3/core/metafunction/iterator.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/core/metafunction/transformation_trait_or.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/view/detail.hpp>
#include <seqan3/range/detail/inherited_iterator_base.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/iterator>
#include <seqan3/std/type_traits>
#include <seqan3/std/view/view_all.hpp>

namespace seqan3::detail
{

// ============================================================================
//  view_take
// ============================================================================

/*!\brief The type returned by seqan3::view::take and seqan3::view::take_or_throw.
 * \tparam urng_t    The type of the underlying ranges, must satisfy seqan3::view::concept.
 * \tparam exactly   Whether to expose sized'ness on the view.
 * \tparam or_throw  Whether to throw an exception when the input is exhausted before the end of line is reached.
 * \implements std::ranges::View
 * \implements std::ranges::RandomAccessRange
 * \implements std::ranges::SizedRange
 * \ingroup view
 *
 * \details
 *
 * Note that most members of this class are generated by ranges::view_interface which is not yet documented here.
 */
template <std::ranges::View urng_t, bool exactly, bool or_throw>
class view_take : public ranges::view_interface<view_take<urng_t, exactly, or_throw>>
{
private:
    //!\brief The underlying range.
    urng_t urange;

    //!\brief The desired target_size.
    size_t target_size;

    //!\brief The sentinel type is identical to that of the underlying range.
    using sentinel_type = sentinel_t<urng_t>;

    //!\brief The iterator type inherits from the underlying type, but overwrites several operators.
    //!\tparam rng_t Should be `urng_t` for defining #iterator and `urng_t const` for defining #const_iterator.
    template <typename rng_t>
    class iterator_type : public inherited_iterator_base<iterator_type<rng_t>, iterator_t<rng_t>>
    {
    private:
        //!\brief The iterator type of the underlying range.
        using base_base_t = iterator_t<rng_t>;
        //!\brief The CRTP wrapper type.
        using base_t      = inherited_iterator_base<iterator_type, iterator_t<rng_t>>;

        //!\brief The current position.
        size_t pos;

        //!\brief The size parameter to the view.
        size_t max_pos;

    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        iterator_type() = default;
        constexpr iterator_type(iterator_type const & rhs) = default;
        constexpr iterator_type(iterator_type && rhs) = default;
        constexpr iterator_type & operator=(iterator_type const & rhs) = default;
        constexpr iterator_type & operator=(iterator_type && rhs) = default;
        ~iterator_type() = default;

        //!\brief Constructor that delegates to the CRTP layer.
        iterator_type(base_base_t const & it) :
            base_t{it}
        {}

        //!\brief Constructor that delegates to the CRTP layer and initialises the members.
        iterator_type(base_base_t it, size_t const _pos, size_t const _max_pos) :
            base_t{it}, pos{_pos}, max_pos(_max_pos)
        {}
        //!\}

        /*!\name Associated types
         * \brief All are derived from the base_base_t.
         * \{
         */
        using difference_type       = typename std::iterator_traits<base_base_t>::difference_type;
        using value_type            = typename std::iterator_traits<base_base_t>::value_type;
        using reference             = typename std::iterator_traits<base_base_t>::reference;
        using pointer               = typename std::iterator_traits<base_base_t>::pointer;
        using iterator_category     = typename std::iterator_traits<base_base_t>::iterator_category;
        //!\}

        /*!\name Arithmetic operators
         * \brief seqan3::detail::inherited_iterator_base operators are used unless specialised here.
         * \{
         */
        iterator_type & operator++() noexcept(noexcept(++base_base_t{}))
        {
            base_base_t::operator++();
            ++pos;
            return *this;
        }

        iterator_type operator++(int) noexcept(noexcept(++base_base_t{}))
        {
            iterator_type cpy{*this};
            ++(*this);
            return cpy;
        }

        iterator_type & operator--() noexcept(noexcept(--base_base_t{}))
        //!\cond
            requires std::BidirectionalIterator<base_base_t>
        //!\endcond
        {
            base_base_t::operator--();
            --pos;
            return *this;
        }

        iterator_type operator--(int) noexcept(noexcept(--base_base_t{}))
        //!\cond
            requires std::BidirectionalIterator<base_base_t>
        //!\endcond
        {
            iterator_type cpy{*this};
            --(*this);
            return cpy;
        }

        iterator_type & operator+=(difference_type const skip) noexcept(noexcept(base_base_t{} += skip))
        //!\cond
            requires std::RandomAccessIterator<base_base_t>
        //!\endcond
        {
            base_base_t::operator+=(skip);
            pos += skip;
            return *this;
        }

        iterator_type & operator-=(difference_type const skip) noexcept(noexcept(base_base_t{} -= skip))
        //!\cond
            requires std::RandomAccessIterator<base_base_t>
        //!\endcond
        {
            base_base_t::operator-=(skip);
            pos -= skip;
            return *this;
        }
        //!\}

        /*!\name Comparison operators
         * \brief We define comparison against self and against the sentinel.
         * \{
         */
        bool operator==(iterator_type const & rhs) const noexcept(!or_throw)
        {
            return static_cast<base_base_t>(*this) == static_cast<base_base_t>(rhs);
        }

        bool operator==(sentinel_type const & rhs) const noexcept(!or_throw)
        {
            if (pos >= max_pos)
                return true;

            if (static_cast<base_base_t>(*this) == rhs)
            {
                if constexpr (or_throw)
                    throw unexpected_end_of_input{"Reached end of input before designated size."};

                return true;
            }
            else
            {
                return false;
            }
        }

        friend bool operator==(sentinel_type const & lhs, iterator_type const & rhs) noexcept(!or_throw)
        {
            return rhs == lhs;
        }

        bool operator!=(sentinel_type const & rhs) const noexcept(!or_throw)
        {
            return !(*this == rhs);
        }

        bool operator!=(iterator_type const & rhs) const noexcept(!or_throw)
        {
            return static_cast<base_base_t>(*this) != static_cast<base_base_t>(rhs);
        }

        friend bool operator!=(sentinel_type const & lhs, iterator_type const & rhs) noexcept(!or_throw)
        {
            return rhs != lhs;
        }
        //!\}

        /*!\name Reference/Dereference operators
         * \brief seqan3::detail::inherited_iterator_base operators are used unless specialised here.
         * \{
         */
        reference operator[](std::make_unsigned_t<difference_type> const n) const noexcept(noexcept(base_base_t{}[0]))
        //!\cond
            requires std::RandomAccessIterator<base_base_t>
        //!\endcond
        {
            pos = n;
            return base_base_t::operator[](n);
        }
        //!\}
    }; // class iterator_type

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The reference_type.
    using reference         = reference_t<urng_t>;
    //!\brief The const_reference type is equal to the reference type if the underlying range is const-iterable.
    using const_reference   = detail::transformation_trait_or_t<seqan3::reference<urng_t const>, void>;
    //!\brief The value_type (which equals the reference_type with any references removed).
    using value_type        = value_type_t<urng_t>;
    //!\brief The size_type is `size_t` if the view is exact, otherwise `void`.
    using size_type         = std::conditional_t<exactly, size_t, void>;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type   = difference_type_t<urng_t>;
    //!\brief The iterator type of this view (a random access iterator).
    using iterator          = iterator_type<urng_t>;
    //!\brief The const_iterator type is equal to the iterator type if the underlying range is const-iterable.
    using const_iterator    = detail::transformation_trait_or_t<std::type_identity<iterator_type<urng_t const>>, void>;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    view_take() = default;
    constexpr view_take(view_take const & rhs) = default;
    constexpr view_take(view_take && rhs) = default;
    constexpr view_take & operator=(view_take const & rhs) = default;
    constexpr view_take & operator=(view_take && rhs) = default;
    ~view_take() = default;

    /*!\brief Construct from another range.
     * \param[in] _urange The underlying range.
     * \param[in] _size   The desired size (after which to stop returning elements).
     * \throws unexpected_end_of_input If `exactly && or_throw && seqan3::sized_range_concept<urng_t>`.
     */
    view_take(urng_t _urange, size_t const _size)
        : urange{std::move(_urange)}, target_size{_size}
    {
        if constexpr (exactly && or_throw && std::ranges::SizedRange<urng_t>)
        {
            if (ranges::size(_urange) < _size)
            {
                throw std::invalid_argument{
                    "You are trying to construct a view::take_exactly_or_throw from a range that is strictly smaller."};
            }
        }
    }
    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the container.
     * \returns Iterator to the first element.
     *
     * If the container is empty, the returned iterator will be equal to end().
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator begin() noexcept
    {
        return {ranges::begin(urange), 0, target_size};
    }

    //!\copydoc begin()
    const_iterator begin() const noexcept
        requires const_iterable_concept<urng_t>
    {
        return {ranges::cbegin(urange), 0, target_size};
    }

    //!\copydoc begin()
    const_iterator cbegin() const noexcept
        requires const_iterable_concept<urng_t>
    {
        return {ranges::cbegin(urange), 0, target_size};
    }

    /*!\brief Returns an iterator to the element following the last element of the range.
     * \returns Iterator to the end.
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    sentinel_type end() noexcept
    {
        return {ranges::end(urange)};
    }

    //!\copydoc end()
    sentinel_type end() const noexcept
        requires const_iterable_concept<urng_t>
    {
        return {ranges::cend(urange)};
    }

    //!\copydoc end()
    sentinel_type cend() const noexcept
        requires const_iterable_concept<urng_t>
    {
        return {ranges::cend(urange)};
    }
    //!\}

    /*!\brief Returns the number of elements in the view (only available for "exactly" specialisation!).
     * \returns The number of elements in the view.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    size_type size() const noexcept
        requires exactly
    {
        return target_size;
    }

    /*!\brief Convert this view into a container implicitly.
     * \tparam container_t Type of the container to convert to; must satisfy seqan3::sequence_container_concept and the
     *                     seqan3::reference_t of both must model std::CommonReference.
     * \returns This view converted to container_t.
     */
    template <sequence_container_concept container_t>
    operator container_t()
    //!\cond
        requires std::CommonReference<reference_t<container_t>, reference>
    //!\endcond
    {
        container_t ret;
        ranges::copy(begin(), end(), std::back_inserter(ret));
        return ret;
    }

    //!\overload
    template <sequence_container_concept container_t>
    operator container_t() const
    //!\cond
        requires std::CommonReference<reference_t<container_t>, reference> && const_iterable_concept<urng_t>
    //!\endcond
    {
        container_t ret;
        ranges::copy(cbegin(), cend(), std::back_inserter(ret));
        return ret;
    }
};

//!\brief Template argument type deduction guide that strips references.
//!\relates seqan3::detail::view_take
template <typename urng_t,
          bool exactly = false,
          bool or_throw = false>
view_take(urng_t, size_t) -> view_take<std::remove_reference_t<urng_t>, exactly, or_throw>;

// ============================================================================
//  take_fn (adaptor definition)
// ============================================================================

/*!\brief View adaptor definition for view::take and view::take_or_throw.
 * \tparam or_throw Whether to throw an exception when the input is exhausted before the end of line is reached.
 */
template <bool exactly, bool or_throw>
class take_fn : public pipable_adaptor_base<take_fn<exactly, or_throw>>
{
private:
    //!\brief Type of the CRTP-base.
    using base_t = pipable_adaptor_base<take_fn<exactly, or_throw>>;

public:
    //!\brief Inherit the base class's Constructors.
    using base_t::base_t;

private:
    //!\brief Befriend the base class so it can call impl().
    friend base_t;

    /*!\brief       Call the view's constructor with the underlying view as argument.
     * \returns     An instance of seqan3::detail::view_take.
     */
    template <std::ranges::View urng_t>
    static auto impl(urng_t urange, size_t const target_size)
    {
        return view_take<urng_t, exactly, or_throw>{std::move(urange), target_size};
    }

    /*!\brief       Call the view's constructor with the underlying range wrapped in seqan3::view::all as argument.
     * \returns     An instance of seqan3::detail::view_take.
     */
    template <std::ranges::ViewableRange urng_t>
    static auto impl(urng_t && urange, size_t const target_size)
    {
        return impl(view::all(std::forward<urng_t>(urange)), target_size);
    }
};

} // namespace seqan3::detail

// ============================================================================
//  view::take (adaptor instance definition)
// ============================================================================

namespace seqan3::view
{

/*!\name General purpose views
 * \{
 */

/*!\brief               A view adaptor that returns the first `size` elements from the underlying range (or less if the
 *                      underlying range is shorter).
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \param[in] size      The target size of the view.
 * \returns             Up to `size` elements of the underlying range.
 * \ingroup view
 *
 * \details
 *
 * ### View properties
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |---------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                                        |
 * | std::ranges::ForwardRange       |                                       | *preserved*                                        |
 * | std::ranges::BidirectionalRange |                                       | *preserved*                                        |
 * | std::ranges::RandomAccessRange  |                                       | *preserved*                                        |
 * | std::ranges::ContiguousRange    |                                       | *preserved*                                        |
 * |                                 |                                       |                                                    |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                       |
 * | std::ranges::View               |                                       | *guaranteed*                                       |
 * | std::ranges::SizedRange         |                                       | *lost*                                             |
 * | std::ranges::CommonRange        |                                       | *lost*                                             |
 * | std::ranges::OutputRange        |                                       | *preserved*                                        |
 * | seqan3::const_iterable_concept  |                                       | *preserved*                                        |
 * |                                 |                                       |                                                    |
 * | seqan3::reference_t             |                                       | seqan3::reference_t<urng_t>                        |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * \snippet test/snippet/range/view/take.cpp usage
 *
 * \hideinitializer
 */
inline auto constexpr take = detail::take_fn<false, false>{};

//!\}

} // namespace seqan3::view
