(ns gibbs
  (:require [clojure.string :as str])
  (:use [clojure.pprint :only [pprint]]))

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
  (vec (keep-indexed
        (fn [i e] (if (not= i n) e))
        v)))

;; (remove-nth [0 1 2 3 4 5] 2)

(defn l-mer [s i len]
  (subs s i (+ i len)))

(defn subseqs [seqs len]
  (vec (for [seq seqs]
         (let [start (rand-int (- (count seq) len))]
           (l-mer seq start len)))))

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

;; ((freq-dist "ATGCC") (nth bases 1))
;; (freq-dist (map #(nth % 0) (subseqs seqs 8)))

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

(defn l-mers
  "Returns all subsequences of s with length len"
  [s len]
  (for [i (range (- (count s) len))]
    (l-mer s i len)))
    
(defn bi [base]
  (.indexOf bases (str base)))

;; (bi \T) => 1

(defn lmer-prob [lmer profile]
  (reduce * (map-indexed (fn [i base]
                           (get-in profile [(bi base) i]))
                         lmer)))

;; from http://rosettacode.org/wiki/Probabilistic_choice#Clojure
(defn to-cdf [pdf]
  (reduce
    (fn [acc n] (conj acc (+ (or (last acc) 0) n)))
    []
    pdf))

; (to-cdf [1 2.0 1 1])

(defn rank-lmers [lmers profile]
  (let [probs      (map #(lmer-prob % profile) lmers)
        nonzeros   (filter #(not= % 0.0) probs)
        max-prob   (apply max probs)
        min-prob   (if (not (zero? (count nonzeros)))
                     (apply min nonzeros))
        prob-dist  (if min-prob
                     (map #(/ % min-prob) probs)
                     (repeat (count lmers) 1))
        prob-cdf   (to-cdf prob-dist)
        rand-f     (rand (last prob-cdf))
        new-start  (count (filter #(< % rand-f) prob-cdf))]
    [new-start max-prob]))

;; this should yield a uniform dist over start positions
(defn test-ranker []
  (rank-lmers ["AAA" "ATG" "CAG"]
              (to-array-2d [[0.5 0.5 0.5]
                            [0.0 0.5 0.0]
                            [0.0 0.0 0.5]
                            [0.5 0.0 0.0]])))

;; it works well
(frequencies
 (for [i (range 10000)]
  (test-ranker)))

(defn gibbs-iter [seqs len ss]
  (let [removed-idx (rand-int (count seqs))
        removed-seq (nth seqs removed-idx)
        new-ss      (remove-nth ss removed-idx)
        profile     (make-profile new-ss len bases)
        [new-start max-prob] (rank-lmers (l-mers removed-seq len)
                                         profile)
        added-ss    (l-mer removed-seq new-start len)]
    ;;(pprint new-ss)
    ;;(pprint profile)
    ;;(pprint max-prob)

    ;; must return new ss in order!
    ;; ss must be a vector for this to work
    [max-prob
     added-ss
     (reduce conj (conj (subvec ss 0 removed-idx)
                        added-ss)
             (subvec ss (inc removed-idx)))]))

;; Idea: stop iterating once we've gone cutoff rounds
;; without seeing an increase in max-prob

(defn gibbs [seqs len cutoff]
  ;; (println "Beginning gibbs")
  ;; rounds = the number of rounds without improvement
  (loop [ss (subseqs seqs len) best 0.0 rounds 0]
    (let [[max-prob added-ss new-ss] (gibbs-iter seqs len ss)
          new-best                   (max max-prob best)
          rounds                     (if (> new-best best)
                                       0
                                       (inc rounds))]
      
      (if (< rounds cutoff)
        (recur new-ss new-best rounds)
        [new-best added-ss new-ss])
      )))

(gibbs seqs 8 50)

(first (sort-by #(* -1 (first %))
                (pmap (fn [_] (gibbs seqs 8 50)) (range 100))))


(defn test-data [n]
  (let [data-str (slurp (str "../data/data" n ".txt"))
        lines    (str/split-lines data-str)]
    (first (sort-by #(* -1 (first %))
                    (pmap (fn [_] (gibbs lines 10 50))
                          (range 20))))))

(test-data 1)
;; [1.0 "TTCGAATTCC" ["TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC"]]
;; [0.00164794921875 "AACAATAT" ["AACAATAT" "ACCTCGCA" "CCGTACTG" "TAAACGAC" "ACACCCTG"]]
;; [1.0 "AATTCGAATT" ["AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT"]]


(test-data 2)
;;[0.19182731573049971 "CTGTCTACTA" ["CTGTCTGCTA" "CTGTCTACTA" "CCGTCTTCTA" "CTGATTACTA" "TCGGCTACTA" "TCAGAGGGAG" "CTGTCAACTA" "CTATCTACTA" "CATATAACAA" "CTGTCTACTA"]]

;; I need to implement relative entropy agains the background
;; should be pretty straightforward.
;; 
