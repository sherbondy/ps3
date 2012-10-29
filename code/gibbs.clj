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

(defn gibbs [seqs len]
  (for [seq seqs]
    (let [start (rand-int (- (count seq) len))
          end   (+ start len)]
      (subs seq start end)))
  )

(gibbs seqs 8)