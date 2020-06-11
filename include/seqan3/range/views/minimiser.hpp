// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Mitra Darvish <mitra.darvish AT fu-berlin.de>
 * \brief Provides seqan3::views::minimiser.
 */

#pragma once

#include <seqan3/std/algorithm>
#include <deque>

#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/type_traits/lazy.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/detail.hpp>

namespace seqan3::detail
{
// ---------------------------------------------------------------------------------------------------------------------
// minimiser_view class
// ---------------------------------------------------------------------------------------------------------------------

/*!\brief The type returned by seqan3::views::minimiser.
 * \tparam urng1_t The type of the underlying ranges, must model std::ranges::forward_range, the reference type must
 *                 model std::totally_ordered. The typical use case is that the reference type is the result of
 *                 seqan3::kmer_hash.
 * \implements std::ranges::view
 * \ingroup views
 *
 * \details
 *
 * Note that most members of this class are generated by ranges::view_interface which is not yet documented here.
 */
template <std::ranges::view urng1_t,
          std::ranges::view urng2_t = std::ranges::empty_view<seqan3::detail::empty_type>>
class minimiser_view : public std::ranges::view_interface<minimiser_view<urng1_t, urng2_t>>
{
private:
    static_assert(std::ranges::forward_range<urng1_t>, "The minimiser_view only works on forward_ranges.");
    static_assert(std::ranges::forward_range<urng2_t>, "The minimiser_view only works on forward_ranges.");
    static_assert(std::totally_ordered<std::ranges::range_reference_t<urng1_t>>,
                  "The reference type of the underlying range must model std::totally_ordered.");

    //!\brief The default argument of the second range.
    using default_urng2_t = std::ranges::empty_view<seqan3::detail::empty_type>;

    //!\brief Boolean variable, which is true, when second range is not of empty type.
    static constexpr bool second_range_is_given = !std::same_as<urng2_t, default_urng2_t>;

    static_assert(!second_range_is_given || std::totally_ordered_with<std::ranges::range_reference_t<urng1_t>,
                                                                      std::ranges::range_reference_t<urng2_t>>,
                  "The reference types of the underlying ranges must model std::totally_ordered_with.");

    //!\brief Whether the given ranges are const_iterable
    static constexpr bool const_iterable = seqan3::const_iterable_range<urng1_t> &&
                                           seqan3::const_iterable_range<urng2_t>;

    //!\brief The first underlying range.
    urng1_t urange1;

    //!\brief The second underlying range.
    urng2_t urange2;

    //!\brief The number of values in one window.
    size_t window_values_size;

    template <typename rng1_t, typename rng2_t>
    class window_iterator;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    minimiser_view() = default; //!< Defaulted.
    minimiser_view(minimiser_view const & rhs) = default; //!< Defaulted.
    minimiser_view(minimiser_view && rhs) = default; //!< Defaulted.
    minimiser_view & operator=(minimiser_view const & rhs) = default; //!< Defaulted.
    minimiser_view & operator=(minimiser_view && rhs) = default; //!< Defaulted.
    ~minimiser_view() = default; //!< Defaulted.

    /*!\brief Construct from a view and a given number of values in one window.
    * \param[in] urange1_ The input range to process. Must model std::ranges::viewable_range and
    *                     std::ranges::forward_range.
    * \param[in] w_ The number of values in one window.
    */
    minimiser_view(urng1_t urange1_, size_t const w_) :
        minimiser_view{std::move(urange1_), default_urng2_t{}, w_}
    {}

    /*!\brief Construct from a non-view that can be view-wrapped and a given number of values in one window.
    * \param[in] urange1_ The input range to process. Must model std::ranges::viewable_range and
    *                     std::ranges::forward_range.
    * \param[in] w_ The number of values in one window.
    */
    template <typename other_urng1_t>
    //!\cond
     requires (std::ranges::viewable_range<other_urng1_t> &&
               std::constructible_from<urng1_t, ranges::ref_view<std::remove_reference_t<other_urng1_t>>>)
    //!\endcond
    minimiser_view(other_urng1_t && urange1_, size_t const w_) :
        urange1{std::views::all(std::forward<other_urng1_t>(urange1_))},
        urange2{default_urng2_t{}},
        window_values_size{w_}
    {}

    /*!\brief Construct from two views and a given number of values in one window.
    * \param[in] urange1_ The first input range to process. Must model std::ranges::viewable_range and
    *                     std::ranges::forward_range.
    * \param[in] urange2_ The second input range to process. Must model std::ranges::viewable_range and
    *                     std::ranges::forward_range.
    * \param[in] w_ The number of values in one window.
    */
    minimiser_view(urng1_t urange1_, urng2_t urange2_, size_t const w_) :
        urange1{std::move(urange1_)},
        urange2{std::move(urange2_)},
        window_values_size{w_}
    {}

    /*!\brief Construct from two non-views that can be view-wrapped and a given number of values in one window.
    * \param[in] urange1_ The input range to process. Must model std::ranges::viewable_range and
    *                     std::ranges::forward_range.
    * \param[in] urange2_ The second input range to process. Must model std::ranges::viewable_range and
    *                     std::ranges::forward_range.
    * \param[in] w_ The number of values in one window.
    */
    template <typename other_urng1_t, typename other_urng2_t>
    //!\cond
     requires (std::ranges::viewable_range<other_urng1_t> &&
               std::constructible_from<urng1_t, ranges::ref_view<std::remove_reference_t<other_urng1_t>>> &&
               std::ranges::viewable_range<other_urng2_t> &&
               std::constructible_from<urng2_t, ranges::ref_view<std::remove_reference_t<other_urng2_t>>>)
    //!\endcond
    minimiser_view(other_urng1_t && urange1_, other_urng2_t && urange2_, size_t const w_) :
        urange1{std::views::all(std::forward<other_urng1_t>(urange1_))},
        urange2{std::views::all(std::forward<other_urng2_t>(urange2_))},
        window_values_size{w_}
    {}
    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the range.
     * \returns Iterator to the first element.
     *
     * \details
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * Strong exception guarantee.
     */
    auto begin()
    {
        return window_iterator<urng1_t, urng2_t>{std::ranges::end(urange1),
                                                 std::ranges::begin(urange1),
                                                 std::ranges::begin(urange2),
                                                 window_values_size};
    }

    //!\copydoc begin()
    auto begin() const
    //!\cond
        requires const_iterable
    //!\endcond
    {
        return window_iterator<urng1_t const, urng2_t const>{std::ranges::end(urange1),
                                                             std::ranges::begin(urange1),
                                                             std::ranges::begin(urange2),
                                                             window_values_size};
    }

    //!\copydoc begin()
    auto cbegin() const
    //!\cond
        requires const_iterable
    //!\endcond
    {
        return begin();
    }

    /*!\brief Returns an iterator to the element following the last element of the range.
     * \returns Iterator to the end.
     *
     * \details
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
    auto end()
    {
        return std::ranges::end(urange1);
    }

    //!\copydoc end()
    auto end() const
    //!\cond
        requires const_iterable
    //!\endcond
    {
        return std::ranges::end(urange1);
    }

    //!\copydoc end()
    auto cend() const
    //!\cond
        requires const_iterable
    //!\endcond
    {
        return end();
    }
    //!\}
};

//!\brief Iterator for calculating minimisers.
template <std::ranges::view urng1_t, std::ranges::view urng2_t>
template <typename rng1_t, typename rng2_t>
class minimiser_view<urng1_t, urng2_t>::window_iterator
{
private:
    //!\brief The sentinel type of the first underlying range.
    using urng1_sentinel_t = std::ranges::sentinel_t<rng1_t>;
    //!\brief The iterator type of the first underlying range.
    using urng1_iterator_t = std::ranges::iterator_t<rng1_t>;
    //!\brief The iterator type of the second underlying range.
    using urng2_iterator_t = std::ranges::iterator_t<rng2_t>;

    template <typename rng1_t_, typename rng2_t_>
    friend class window_iterator;

public:
    /*!\name Associated types
     * \{
     */
    //!\brief Type for distances between iterators.
    using difference_type = std::ranges::range_difference_t<rng1_t>;
    //!\brief Value type of this iterator.
    using value_type = std::ranges::range_value_t<rng1_t>;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief Reference to `value_type`.
    using reference = value_type;
    //!\brief Tag this class as a forward iterator.
    using iterator_category = std::forward_iterator_tag;
    //!\brief Tag this class as a forward iterator.
    using iterator_concept = iterator_category;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    window_iterator() = default; //!< Defaulted.
    window_iterator(window_iterator const &) = default; //!< Defaulted.
    window_iterator(window_iterator &&) = default; //!< Defaulted.
    window_iterator & operator=(window_iterator const &) = default; //!< Defaulted.
    window_iterator & operator=(window_iterator &&) = default; //!< Defaulted.
    ~window_iterator() = default; //!< Defaulted.

    //!\brief Allow iterator on a const range to be constructible from an iterator over a non-const range.
    template <typename rng1_t_, typename rng2_t_>
    //!\cond
        requires (std::same_as<std::remove_const_t<urng1_t>, rng1_t_> &&
                  std::same_as<std::remove_const_t<urng2_t>, rng2_t_>)
     //!\endcond
    window_iterator(window_iterator<rng1_t_, rng2_t_> it) :
        minimiser_value{std::move(it.minimiser_value)},
        urng1_sentinel{std::move(it.urng1_sentinel)},
        urng1_iterator{std::move(it.urng1_iterator)},
        urng2_iterator{std::move(it.urng2_iterator)},
        window_values{std::move(it.window_values)}
    { }

    /*!\brief Construct from begin and end iterators of a given range over std::totally_ordered values, and the number
              of values per window.
    * \param[in] urng1_iterator Iterator pointing to the first position of the first std::totally_ordered range.
    * \param[in] urng1_sentinel Iterator pointing to the last position of the first std::totally_ordered range.
    * \param[in] urng2_iterator Iterator pointing to the first position of the second std::totally_ordered range.
    * \param[in] window_values_size The number of values in one window.
    *
    * \details
    *
    * Looks at the number of values per window in two ranges, returns the smallest between both as minimiser and
    * shifts then by one to repeat this action. If a minimiser in consecutive windows is the same, it is returned only
    * once.
    *
    */
    window_iterator(urng1_sentinel_t urng1_sentinel,
                    urng1_iterator_t urng1_iterator,
                    urng2_iterator_t urng2_iterator,
                    size_t window_values_size) :
        urng1_sentinel{std::move(urng1_sentinel)},
        urng1_iterator{std::move(urng1_iterator)},
        urng2_iterator{std::move(urng2_iterator)}
    {
        size_t size = std::ranges::distance(urng1_iterator, urng1_sentinel);
        window_values_size = std::min<size_t>(window_values_size, size);

        window_first(window_values_size);
    }
    //!\}

    //!\anchor window_iterator_comparison
    //!\name Comparison operators
    //!\{

    //!\brief Compare to another window_iterator.
    friend bool operator==(window_iterator const & lhs, window_iterator const & rhs)
    {
        return (lhs.urng1_iterator == rhs.urng1_iterator) &&
               (lhs.window_values.size() == rhs.window_values.size());
    }

    //!\brief Compare to another window_iterator.
    friend bool operator!=(window_iterator const & lhs, window_iterator const & rhs)
    {
        return !(lhs == rhs);
    }

    //!\brief Compare to the sentinel of the underlying range.
    friend bool operator==(window_iterator const & lhs, urng1_sentinel_t const &)
    {
        return lhs.urng1_iterator == lhs.urng1_sentinel;
    }

    //!\brief Compare to the sentinel of the underlying range.
    friend bool operator==(urng1_sentinel_t const & lhs, window_iterator const & rhs)
    {
        return rhs == lhs;
    }

    //!\brief Compare to the sentinel of the underlying range.
    friend bool operator!=(urng1_sentinel_t const & lhs, window_iterator const & rhs)
    {
        return !(lhs == rhs);
    }

    //!\brief Compare to the sentinel of the underlying range.
    friend bool operator!=(window_iterator const & lhs, urng1_sentinel_t const & rhs)
    {
        return !(lhs == rhs);
    }
    //!\}

    //!\brief Pre-increment.
    window_iterator & operator++() noexcept
    {
        next_unique_minimiser();
        return *this;
    }

    //!\brief Post-increment.
    window_iterator operator++(int) noexcept
    {
        window_iterator tmp{*this};
        next_unique_minimiser();
        return tmp;
    }

    //!\brief Return the minimiser.
    value_type operator*() const noexcept
    {
        return minimiser_value;
    }

private:
    //!\brief The minimiser value.
    value_type minimiser_value{};

    //!brief Iterator to last element in range.
    urng1_sentinel_t urng1_sentinel;

    //!\brief Iterator to the rightmost value of one window.
    urng1_iterator_t urng1_iterator;

    //!\brief Iterator to the rightmost value of one window of the second range.
    urng2_iterator_t urng2_iterator;

    //!\brief Stored values per window. It is necessary to store them, because a shift can remove the current minimiser.
    std::deque<value_type> window_values{};

    //!\brief Increments iterator by 1.
    void next_unique_minimiser()
    {
        while (!next_minimiser()) {}
    }

    //!\brief Returns new window value, when only one range is given.
    auto window_value()
    //!\cond
        requires (!second_range_is_given)
    //!\endcond
    {
        return *urng1_iterator;
    }

    //!\brief Returns new window value, when two ranges are given.
    auto window_value()
    //!\cond
        requires second_range_is_given
    //!\endcond
    {
        return std::min(*urng1_iterator, *urng2_iterator);
    }

    //!\brief Advance first range.
    void advance_window()
    //!\cond
        requires (!second_range_is_given)
    //!\endcond
    {
        ++urng1_iterator;
    }

    //!\brief Advance both ranges.
    void advance_window()
    //!\cond
        requires second_range_is_given
    //!\endcond
    {
        ++urng1_iterator;
        ++urng2_iterator;
    }

    //!\brief Calculates minimisers for the first window.
    void window_first(size_t const window_values_size)
    {
        if (window_values_size == 0u)
            return;

        for (size_t i = 0u; i < window_values_size - 1u; ++i)
        {
            window_values.push_back(window_value());
            advance_window();
        }
        window_values.push_back(window_value());
        minimiser_value = *std::ranges::min_element(window_values);
    }

    /*!\brief Calculates the next minimiser value.
     *
     * \details
     * For the following windows, we remove the first window value (is now not in window_values) and add the new
     * value that results from the window shifting.
     */
    bool next_minimiser()
    {
        advance_window();
        if (urng1_iterator == urng1_sentinel)
            return true;

        value_type const new_value = window_value();
        value_type const first_window_value = window_values.front();

        window_values.pop_front();
        window_values.push_back(new_value);

        if (minimiser_value == first_window_value)
        {
            minimiser_value = *std::ranges::min_element(window_values);
            return true;
        }

        if (new_value < minimiser_value)
        {
            minimiser_value = new_value;
            return true;
        }

        return false;
    }
};

//!\brief A deduction guide for the view class template.
template <std::ranges::viewable_range rng1_t>
minimiser_view(rng1_t &&, size_t const window_values_size) -> minimiser_view<std::views::all_t<rng1_t>>;

//!\brief A deduction guide for the view class template.
template <std::ranges::viewable_range rng1_t, std::ranges::viewable_range rng2_t>
minimiser_view(rng1_t &&, rng2_t &&, size_t const window_values_size) -> minimiser_view<std::views::all_t<rng1_t>,
                                                                                        std::views::all_t<rng2_t>>;

// ---------------------------------------------------------------------------------------------------------------------
// minimiser_fn (adaptor definition)
// ---------------------------------------------------------------------------------------------------------------------

//![adaptor_def]
//!\brief views::minimiser's range adaptor object type (non-closure).
struct minimiser_fn
{
    //!\brief Store the number of values in one window and return a range adaptor closure object.
    constexpr auto operator()(size_t const window_values_size) const
    {
        return adaptor_from_functor{*this, window_values_size};
    }

    /*!\brief Call the view's constructor with two arguments: the underlying view and an integer indicating how many
     *        values one window contains.
     * \param[in] urange1 The input range to process. Must model std::ranges::viewable_range and
     *                    std::ranges::forward_range.
     * \param[in] window_values_size The number of values in one window.
     * \returns  A range of converted values.
     */
    template <std::ranges::range urng1_t>
    constexpr auto operator()(urng1_t && urange1, size_t const window_values_size) const
    {
        static_assert(std::ranges::viewable_range<urng1_t>,
                      "The range parameter to views::minimiser cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng1_t>,
                      "The range parameter to views::minimiser must model std::ranges::forward_range.");

        if (window_values_size == 1) // Would just return urange1 without any changes
            throw std::invalid_argument{"The chosen window_values_size is not valid. "
                                        "Please choose a value greater than 1 or use two ranges."};

        return minimiser_view{urange1, window_values_size};
    }


    //!\brief Store the number of values in one window and the second range and return a range adaptor closure object.
    template <std::ranges::range urng2_t>
    constexpr auto operator()(size_t const window_values_size, urng2_t && urange2) const
    {
        return adaptor_from_functor{*this, window_values_size, urange2};
    }

    /*!\brief Call the view's constructor with two arguments: the underlying view and an integer indicating how many
     *        values one window contains.
     * \param[in] urange1 The input range to process. Must model std::ranges::viewable_range and
     *                    std::ranges::forward_range.
     * \param[in] window_values_size The number of values in one window.
     * \param[in] urange2 The second input range to process. Must model std::ranges::viewable_range and
     *                    std::ranges::forward_range.
     * \returns A range of converted values.
     */
    template <std::ranges::range urng1_t, std::ranges::range urng2_t>
    constexpr auto operator()(urng1_t && urange1, size_t const window_values_size, urng2_t && urange2) const
    {
        static_assert(std::ranges::viewable_range<urng1_t>,
                      "The range1 parameter to views::minimiser cannot be a temporary of a non-view range.");
        static_assert(std::ranges::viewable_range<urng2_t>,
                      "The range2 parameter to views::minimiser cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng1_t>,
                      "The range1 parameter to views::minimiser must model std::ranges::forward_range.");
        static_assert(std::ranges::forward_range<urng2_t>,
                      "The range2 parameter to views::minimiser must model std::ranges::forward_range.");

        if (std::ranges::size(urange1) != std::ranges::size(urange2))
            throw std::invalid_argument{"The two ranges do not have the same size."};

        return minimiser_view{urange1, urange2, window_values_size};
    }
};
//![adaptor_def]

} // namespace seqan3::detail

namespace seqan3::views
{

/*!\name General purpose views
 * \{
 */

/*!\brief Computes minimisers for a range of comparable values. A minimiser is the smallest value in a window.
 * \tparam urng1_t The type of the first range being processed. See below for requirements. [template
 *                 parameter is omitted in pipe notation]
 * \tparam urng2_t The type of the range second being processed. See below for requirements. [template
 *                 parameter is omitted in pipe notation]
 * \param[in] urange1 The range being processed. [parameter is omitted in pipe notation]
 * \param[in] window_values_size The number of values in one window.
 * \returns A range of std::totally_ordered where each value is the minimal value for one window. See below for the
 *          properties of the returned range.
 * \ingroup views
 *
 * \details
 *
 * A minimiser is the smallest value in a window. For example for the following list of hash values
 * [28, 100, 9, 23, 4, 1, 72, 37, 8] and 4 as `window_values_size`, the minimiser values are [9,4,1]. The minimiser can
 * be calculated for one given range or for two given ranges, where the minimizer is the smallest value in both windows.
 * For example for the following list of hash values [28, 100, 9, 23, 4, 1, 72, 37, 8] and
 * [30, 2, 11, 101, 199, 73, 34, 900] and 4 as `window_values_size`, the minimiser values are [2,4,1].
 *
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng1_t` (underlying range type)  | `rrng_t` (returned range type)   |
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
 * | seqan3::const_iterable_range     |                                    | *preserved*¹                     |
 * |                                  |                                    |                                  |
 * | std::ranges::range_reference_t   | std::totally_ordered               | std::totally_ordered             |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ¹ The marked properties are only *preserved* if only one range is given or both ranges are of type
 * seqan3::const_iterable_range.
 *
 * ### Example
 *
 * \include test/snippet/range/views/minimiser.cpp
 *
 * \hideinitializer
 */
inline constexpr auto minimiser = detail::minimiser_fn{};

//!\}

} // namespace seqan3::views
