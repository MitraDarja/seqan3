// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides various auxiliary functions with which parts of the configurations can be checked.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <array>
#include <type_traits>

#include <seqan3/core/configuration/detail/concept.hpp>
#include <seqan3/utility/tuple/concept.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// Type trait is_configuration_valid
// ----------------------------------------------------------------------------

/*!\brief Checks if a given type is compatible with a list of other types.
 * \implements seqan3::unary_type_trait
 * \ingroup algorithm
 * \tparam query_t       The type to check for compatibility.
 * \tparam compare_types The types to compare against.
 *
 * \details
 *
 * Checks if the type is from the same algorithm configuration and if it can be combined with any of the
 * existing elements in the current configuration.
 */
template <config_element query_t, config_element ... compare_types>
struct is_configuration_valid :
    public std::conditional_t<
        (std::is_same_v<std::remove_cvref_t<decltype(query_t::id)>, std::remove_cvref_t<decltype(compare_types::id)>> && ...) &&
        (compatibility_table<std::remove_cvref_t<decltype(query_t::id)>>
                [static_cast<std::underlying_type_t<std::remove_cvref_t<decltype(query_t::id)>>>(query_t::id)]
                [static_cast<std::underlying_type_t<std::remove_cvref_t<decltype(query_t::id)>>>(compare_types::id)] && ...),
        std::true_type,  // If condition is true.
        std::false_type  // If condition is false.
    >
{};

/*!\brief Helper variable template to check for valid configuration composites (unary_type_trait shortcut).
 * \relates seqan3::detail::is_configuration_valid
 * \ingroup algorithm
 */
template <typename query_t, typename ...compare_types>
inline constexpr bool is_configuration_valid_v = is_configuration_valid<query_t, compare_types...>::value;

// ----------------------------------------------------------------------------
// Metafunction is_same_configuration_f
// ----------------------------------------------------------------------------

/*!\brief Helper meta function to check if a template type is contained in a seqan3::configuration.
 * \ingroup algorithm
 *
 * \details
 *
 * This helper meta function is used to provide the `get` and `get_or` interface for template template types.
 */
template <template <typename ...> typename query_t>
struct is_same_configuration_f
{
    /*!\brief A type template that evaluates to std::true_type if the given type is a specialization of `query_t`,
     *        otherwise std::false_type.
     * \tparam compare_type The type to compare against `query_t`.
     */
    template <typename compare_type>
    using invoke = is_type_specialisation_of<compare_type, query_t>;
};

} // namespace seqan3::detail
