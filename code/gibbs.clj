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

(defn log2 [n]
  (/ (Math/log n) (Math/log 2)))

(defn rel-entropy [len profile bkgd]
  (reduce + (for [j (range len)]
              (reduce + (for [r bases]
                          (let [p-rj (get-in profile [(bi r) j])]
                            (if (> p-rj 0.0)
                              (* p-rj (log2 (/ p-rj (bkgd r))))
                              0.0)))))))

;; from http://rosettacode.org/wiki/Probabilistic_choice#Clojure
(defn to-cdf [pdf]
  (reduce
    (fn [acc n] (conj acc (+ (or (last acc) 0) n)))
    []
    pdf))

; (to-cdf [1 2.0 1 1])



(defn rank-lmers [lmers profile & [bkgd]]
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

(defn gibbs-iter [seqs len ss & [bkgd]]
  (let [removed-idx (rand-int (count seqs))
        removed-seq (nth seqs removed-idx)
        new-ss      (remove-nth ss removed-idx)
        profile     (make-profile new-ss len bases)
        entropy     (if bkgd
                      (rel-entropy len profile bkgd)
                      0.0)
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
             (subvec ss (inc removed-idx)))
     entropy]))

;; Idea: stop iterating once we've gone cutoff rounds
;; without seeing an increase in max-prob

(defn gibbs [seqs len cutoff & [bkgd]]
  ;; (println "Beginning gibbs")
  ;; rounds = the number of rounds without improvement
  (loop [ss (subseqs seqs len) best 0.0 rounds 0]
    (let [[max-prob added-ss new-ss entropy] (gibbs-iter seqs len ss bkgd)
          metric                     (if bkgd entropy max-prob)
          new-best                   (max metric best)
          rounds                     (if (> new-best best)
                                       0
                                       (inc rounds))]
      
      (if (< rounds cutoff)
        (recur new-ss new-best rounds)
        [max-prob added-ss new-ss entropy])
      )))

;; (gibbs seqs 8 50)

;; (first (sort-by #(* -1 (first %))
;;                (pmap (fn [_] (gibbs seqs 8 50)) (range 100))))


(defn test-data [n & [bkgd]]
  (let [data-str (slurp (str "../data/data" n ".txt"))
        lines    (str/split-lines data-str)
        sortfn   (if bkgd last first)]
    (sort-by #(* -1 (sortfn %))
             (pmap (fn [_] (gibbs lines 10 50 bkgd))
                   (range 20)))))

(test-data 1)
;; [1.0 "TTCGAATTCC" ["TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC"]]
;; [1.0 "AATTCGAATT" ["AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT"]]

(def uniform-bases {\A 0.25 \T 0.25 \C 0.25 \G 0.25})
(test-data 1 uniform-bases)
;; ([1.0 "ATTCGAATTC" ["ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC"] 20.0]

(test-data 2)
;;[0.19182731573049971 "CTGTCTACTA" ["CTGTCTGCTA" "CTGTCTACTA" "CCGTCTTCTA" "CTGATTACTA" "TCGGCTACTA" "TCAGAGGGAG" "CTGTCAACTA" "CTATCTACTA" "CATATAACAA" "CTGTCTACTA"]]
;; [0.015109856187613529 "CGTCTACCAG" ["TGTCTACTAA" "TGTCTACTAA" "CCCCCCAAGC" "TTTCTTTCAA" "TGTCAACTAG" "CGTCTTATAA" "CGTCTACCAG" "TTTCTTCCTA" "TGTCAACTAA" "CTTCATACAA"] 0.0]

(test-data 2 uniform-bases)
;; [0.09293977287302675 "CAACGGCTCT" ["TCTGTCTGCT" "TCGGTGTACT" "CAACGGCTCT" "TCGGTCTACC" "TAAGTCTACT" "TCAGTCTACT" "TCTGTCTGCT" "TCTATCTACT" "TCGGTCTACT" "TCTGTCTACT"] 15.692269256226584]
;; [0.007929370626154036 "TTTGTCTACT" ["TTTGTCTACT" "GCTGTCTACT" "GTCTTCTACT" "AATATCAACT" "TAAGTCTACT" "TCAGTCTACT" "TCTGTCAACT" "AAAAGGACGT" "TCTGTCAACT" "TCTGTCTACT"] 11.632687891823114]

;; By supplying bkgd, we try to maximize relative entropy during
;; the final stage

(test-data 3)
;; [0.11954769052989207 "AAAAAAAAAA" ["AAAAAAAGAA" "AAAAAAAAGA" "AAAGGATAGA" "AAAAATAGAA" "AAAAAAAAAG" "TAAGAATAAA" "AAAAAATTAA" "AAAAAAAAAA" "TAAAGATAAA" "TAACCAGTGG" "AAAAAATAAA" "TAAAAAAAAA" "AAAGAAAAGA" "AAAAGAAAGA" "AAAAGTAAGA" "AAAGATTATA" "AAAAATTTAA" "AAAAAATAAA" "AAAAAAAAAA" "TAAGAACAAA" "TAAAATATAA" "TATAAAAGGA" "AAAAAAAGAA" "AAAGAATAAA" "AAAAAAAAAA" "TAAAATAATA" "AAAGAAAAAA" "TAAAAATATA" "TAAAGAAGGA" "TATAAATATA" "AAAAGAAAAA" "AAAAAACAAA" "AAAAAAAGAA" "AAAAAAAAGA"] 0.0]

;; No good!

(def low-gc-dist {\A 0.315 \T 0.315 \C 0.185 \G 0.185})
(test-data 3 low-gc-dist)
;; [0.013284171165914994 "AAGAATAAAA" ["AAAAAGAAAA" "AAAAAAAAAA" "AGGAAGAAAT" "AAAAAATAAA" "AAAAAGAAAA" "AAGAATAAAA" "AAATTAAAAA" "AAAAAAAAAA" "AAGATAAAAA" "CCACAGAATA" "AAAAATAAAA" "AAATAATAAT" "AGGTAGAAAT" "AGAATATAAT" "AGAATAAAAA" "AAGTAGAAGT" "AAAATGAAAA" "AAAAATAAAA" "AAAAAAAAAA" "AAGAATAAAA" "AAAATATAAT" "AAAAAAAAAA" "AAGAAAAAAA" "AAGAATAAAA" "AAAAAAAAAA" "AAAAAGAAGA" "AAAAAGAAAA" "AAGAAGAAAA" "AAATTGAAAA" "AAATATAAAT" "AAAAAGAAAA" "AAAAAAAAAA" "AAGAAGAAAA" "AAAAAAAAAA"] 10.745940720821562]


(test-data 4 low-gc-dist)
;; [1.2642597750635616E-4 "TTGCTTTTAT" ["TTTTTTTTTT" "TTTTTGTTTT" "TTTATTGTTA" "TTTTTTTTTT" "TTTTTTTTTT" "TTTCCATTTT" "TTTTTTTTTA" "TTTCTTTTTG" "TATCTAATTT" "TATTTTTTTT" "TTTTTTTTTT" "TTTCCTTCTT" "TTTACTATTT" "TTTTTTTCAT" "TTTTTTGTTG" "TTTACTGTAT" "TTTTTTTTAT" "TTTTTTGTTA" "TTTTTTTTTT" "TTTTTTGTTT" "TTTTTTTTTT" "TTTCCTGTAT" "TTTATTATTT" "TTTCCTTTTT" "TTTTCTTTTT" "TTTTTTTCTT" "TTTCTTCTTT" "TTTTTTTCTT" "TTTTCTGTTT" "TTTTTTTTTT" "TCTATTTTTT" "TTTTTTATTT" "TTTTCTTTTT" "TTTTTATTTT" "TTGTTTATTG" "TTTTCATTAT" "TTTATGCTTT" "TTTTTTGTTT" "TTTTTTTTTT" "TTTTCTTCTT" "TTTCCTCTTT" "TTTTCGGTTT" "TTTACTTTTG" "TTTTCATTTT" "TTTTTTTTTT" "TTTACTTCTT" "TCTATTGTTT" "TTTACTTTTT" "TTTTCAATTT" "TTTACTTTTT" "TTGCTTTTAT" "TTGTTTTTTG" "TTTCCTTTTT" "TTTTCATTTT" "TTTTTTTTTT" "TTTTTTTTTT"] 10.634170516353091]
