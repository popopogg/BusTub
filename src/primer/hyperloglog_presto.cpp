//===----------------------------------------------------------------------===//
//
//                         BusTub
//
// hyperloglog_presto.cpp
//
// Identification: src/primer/hyperloglog_presto.cpp
//
// Copyright (c) 2015-2025, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#include "primer/hyperloglog_presto.h"
#include <bitset>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <shared_mutex>
#include "common/util/hash_util.h"
#include "primer/hyperloglog.h"

namespace bustub {

/** @brief Parameterized constructor. */
template <typename KeyType>
HyperLogLogPresto<KeyType>::HyperLogLogPresto(int16_t n_leading_bits) : cardinality_(0) {
  if (n_leading_bits <= 0) {
    b_ = 0;
  } else {
    b_ = n_leading_bits;
  }
  m_ = (1 << b_);
  dense_bucket_.resize(m_);
}

/** @brief Element is added for HLL calculation. */
template <typename KeyType>
auto HyperLogLogPresto<KeyType>::AddElem(KeyType val) -> void {
  /** @TODO(student) Implement this function! */
  hash_t hs = CalculateHash(val);
  std::bitset<BITSET_CAPACITY> binary(hs);
  size_t j = (binary >> (BITSET_CAPACITY - b_)).to_ullong();
  int64_t tot = 0;
  // 查找从右到左的0的数量
  auto count_trailing_zeros = [this](const std::bitset<64> &bits) {
    size_t count = 0;
    const size_t max_check = bits.size() - b_;
    // 保持原逻辑的循环边界，避免无符号整数下溢
    for (size_t i = 0; i < max_check; ++i) {
      if (!bits[i]) {  // 当前位为0则计数
        count++;
      } else {  // 遇到1则终止
        break;
      }
    }
    // 处理全0的特殊情况（保持原逻辑）
    if (count == max_check) {
      count = max_check;
    }
    return count;
  };
  tot = count_trailing_zeros(binary);
  int64_t old_value = dense_bucket_[j].to_ullong();
  // if overflow
  if (overflow_bucket_.find(j) != overflow_bucket_.end()) {
    old_value += (overflow_bucket_[j].to_ulong()) << DENSE_BUCKET_SIZE;
  }
  int64_t new_value = std::max(old_value, tot);
  auto overflow_val = new_value >> DENSE_BUCKET_SIZE;
  // write lock
  std::unique_lock<std::shared_mutex> lk(rwmtx_);

  if (overflow_val > 0) {
    overflow_bucket_[j] = overflow_val;
    dense_bucket_[j] = new_value - (overflow_val << DENSE_BUCKET_SIZE);
    return;
  }
  dense_bucket_[j] = new_value;
}

/** @brief Function to compute cardinality. */
template <typename T>
auto HyperLogLogPresto<T>::ComputeCardinality() -> void {
  /** @TODO(student) Implement this function! */
  // read lock
  std::shared_lock<std::shared_mutex> lk(rwmtx_);

  double sum = 0.0;
  int m = dense_bucket_.size();
  if (m == 0) {
    return;
  }
  for (int j = 0; j < m; ++j) {
    int64_t val = dense_bucket_[j].to_ullong();
    if (overflow_bucket_.find(j) != overflow_bucket_.end()) {
      val += overflow_bucket_[j].to_ullong() << DENSE_BUCKET_SIZE;
    }
    sum += 1.0 / std::pow(2, val);
  }
  cardinality_ = static_cast<size_t>(std::floor(CONSTANT * m * m / sum));
}

template class HyperLogLogPresto<int64_t>;
template class HyperLogLogPresto<std::string>;
}  // namespace bustub
