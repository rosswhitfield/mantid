// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidKernel/StringTokenizer.h"
#include "MantidKernel/Strings.h"
#include <algorithm>
#include <iterator> //cbegin,cend
#include <stdexcept>

namespace {

// If the final character is a separator, we need to add an empty string to
// tokens.
void addEmptyFinalToken(const std::string &str, const std::string &delims, std::vector<std::string> &tokens) {

  const auto pos = std::find(delims.cbegin(), delims.cend(), str.back());

  if (pos != delims.cend()) {
    tokens.emplace_back();
  }
}

// generic tokenizer using std::find_first_of modelled after
// http://tcbrindle.github.io/posts/a-quicker-study-on-tokenising/
// MIT licensed.
template <class InputIt, class ForwardIt, class BinOp>
void for_each_token(InputIt first, InputIt last, ForwardIt s_first, ForwardIt s_last, BinOp binary_op) {
  while (first != last) {
    const auto pos = std::find_first_of(first, last, s_first, s_last);
    binary_op(first, pos);
    if (pos == last)
      break;
    first = std::next(pos);
  }
}

void splitKeepingWhitespaceEmptyTokens(const std::string &str, const std::string &delims,
                                       std::vector<std::string> &output) {
  output.clear();
  for_each_token(str.cbegin(), str.cend(), delims.cbegin(), delims.cend(),
                 [&output](std::string::const_iterator first, std::string::const_iterator second) {
                   output.emplace_back(first, second);
                 });
}

void splitKeepingWhitespaceIgnoringEmptyTokens(const std::string &str, const std::string &delims,
                                               std::vector<std::string> &output) {
  output.clear();
  for_each_token(str.cbegin(), str.cend(), delims.cbegin(), delims.cend(),
                 [&output](std::string::const_iterator first, std::string::const_iterator second) {
                   if (first != second)
                     output.emplace_back(first, second);
                 });
}

void splitIgnoringWhitespaceKeepingEmptyTokens(const std::string &str, const std::string &delims,
                                               std::vector<std::string> &output) {
  output.clear();
  for_each_token(str.cbegin(), str.cend(), delims.cbegin(), delims.cend(),
                 [&output](std::string::const_iterator first, std::string::const_iterator second) {
                   output.emplace_back(first, second);
                   Mantid::Kernel::Strings::stripInPlace(output.back());
                 });
}

void splitIgnoringWhitespaceEmptyTokens(const std::string &str, const std::string &delims,
                                        std::vector<std::string> &output) {
  output.clear();
  for_each_token(str.cbegin(), str.cend(), delims.cbegin(), delims.cend(),
                 [&output](std::string::const_iterator first, std::string::const_iterator second) {
                   if (first != second) {
                     output.emplace_back(first, second);
                     Mantid::Kernel::Strings::stripInPlace(output.back());
                     if (output.back().empty())
                       output.pop_back();
                   }
                 });
}
} // namespace

/**
 * Constructor requiring a string to tokenize and a string of separators.
 * @param str Input string to be separated into tokens.
 * @param separators List of characters used to separate the input string.
 * @param options  tokenizer settings. The number can be found using the
 * StringTokenizer::Options enum
 * @throw Throws std::runtime_error if options > 7.
 */
Mantid::Kernel::StringTokenizer::StringTokenizer(const std::string &str, const std::string &separators,
                                                 unsigned options) {

  // if str is empty, then there is no work to do. exit early.
  if (str.empty())
    return;

  // see comments above for the different options split0,split1,split2 and
  // split3 implement.
  // cases 0-3 will check for a separator in the last place and insert an empty
  // token at the end.
  // cases 4-7 will not check and ignore a potential empty token at the end.
  switch (options) {
  case 0:
    splitKeepingWhitespaceEmptyTokens(str, separators, m_tokens);
    addEmptyFinalToken(str, separators, m_tokens);
    return;
  case TOK_IGNORE_EMPTY:
    splitKeepingWhitespaceIgnoringEmptyTokens(str, separators, m_tokens);
    return;
  case TOK_TRIM:
    splitIgnoringWhitespaceKeepingEmptyTokens(str, separators, m_tokens);
    addEmptyFinalToken(str, separators, m_tokens);
    return;
  case (TOK_TRIM | TOK_IGNORE_EMPTY):
    splitIgnoringWhitespaceEmptyTokens(str, separators, m_tokens);
    return;
  case TOK_IGNORE_FINAL_EMPTY_TOKEN:
    splitKeepingWhitespaceEmptyTokens(str, separators, m_tokens);
    return;
  case (TOK_IGNORE_FINAL_EMPTY_TOKEN | TOK_IGNORE_EMPTY):
    splitKeepingWhitespaceIgnoringEmptyTokens(str, separators, m_tokens);
    return;
  case (TOK_IGNORE_FINAL_EMPTY_TOKEN | TOK_TRIM):
    splitIgnoringWhitespaceKeepingEmptyTokens(str, separators, m_tokens);
    return;
  case (TOK_IGNORE_FINAL_EMPTY_TOKEN | TOK_TRIM | TOK_IGNORE_EMPTY):
    splitIgnoringWhitespaceEmptyTokens(str, separators, m_tokens);
    return;
  }

  // This point is reached only if options > 7.
  throw std::runtime_error("Invalid option passed to Mantid::Kernel::StringTokenizer:" + std::to_string(options));
}
