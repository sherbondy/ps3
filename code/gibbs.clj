(ns gibbs
  (:require [clojure.string :as str]))

;; Much thanks to Neil C. Jones and Pavel A. Pevzner
;; for writing An Introduction to Bioinformatics Algorithms

;; This is a Gibbs sampler for motif discovery

(def seqs
  ["GTAAACAATATTTATAGC"
   "AAAATTTACCTCGCAAGG"
   "CCGTACTGTCAAGCGTGG"
   "TGAGTAAACGACGTCCCA"
   "TACTTAACACCCTGTCAA"])

(defn remove-nth [v n]
  (keep-indexed
   (fn [i e] (if (not= i n) e))
   v))

(defn subseqs [seqs len]
  (for [seq seqs]
    (let [start (rand-int (- (count seq) len))
          end   (+ start len)]
      (subs seq start end))))  

(defn gibbs [seqs len]
  (let [ss (subseqs seqs len)
        sample-idx (rand-int (count seqs))
        new-ss (remove-nth ss sample-idx)]
    [sample-idx new-ss]))

(gibbs seqs 8)