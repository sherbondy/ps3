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

(def bases "ATGC")

(defn remove-nth [v n]
  (keep-indexed
   (fn [i e] (if (not= i n) e))
   v))

(defn subseqs [seqs len]
  (for [seq seqs]
    (let [start (rand-int (- (count seq) len))
          end   (+ start len)]
      (subs seq start end))))

;; (str (get-in seqs [0 0]))

(defn freq-dist
  "Returns the proportionate frequency distribution
   over the elements in a sequence. E.g.
   (freq-dist \"ATGCC\") returns
   {A 0.2, T 0.2, G 0.2, C 0.4}"
  [v]
  (let [fcount (float (count v))]
    (into {}
          (map
           (fn [[k v]] {k (/ v fcount)})
           (frequencies v)))))

((freq-dist "ATGCC") (nth bases 1))

(freq-dist (map #(nth % 0) (subseqs seqs 8)))

(defn make-profile [subseqs size bases]
  ;; first index = letter
  (let [arr (make-array Float/TYPE 4 size)
        bcount (count bases)]
    (doseq [i (range size)]
      (let [freqs (freq-dist (map #(nth % i) subseqs))]
        (doseq [j (range bcount)]
          (aset-float arr j i
                      (get freqs (nth bases j) 0.0))
      )))
    arr))

(defn gibbs [seqs len]
  (let [ss (subseqs seqs len)
        sample-idx (rand-int (count seqs))
        new-ss (remove-nth ss sample-idx)
        profile (make-profile new-ss len bases)]
    [sample-idx new-ss (aget profile 0 0)]))

(gibbs seqs 8)