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
     entropy
     added-ss
     (reduce conj (conj (subvec ss 0 removed-idx)
                        added-ss)
             (subvec ss (inc removed-idx)))]))

;; Idea: stop iterating once we've gone cutoff rounds
;; without seeing an increase in max-prob

(defn gibbs [seqs len cutoff & [bkgd]]
  ;; (println "Beginning gibbs")
  ;; rounds = the number of rounds without improvement
  (loop [ss (subseqs seqs len) best 0.0 rounds 0 best-ss []]
    (let [[max-prob entropy added-ss new-ss] (gibbs-iter seqs len ss bkgd)
          metric                     (if bkgd entropy max-prob)
          new-best                   (max metric best)
          best-ss                    (if (= new-best metric) new-ss best-ss)
          rounds                     (if (> new-best best)
                                       0
                                       (inc rounds))]
      
      (if (< rounds cutoff)
        (recur new-ss new-best rounds best-ss)
        {:best-p new-best :best-ss best-ss})
      )))

;; (gibbs seqs 8 50)
;; {:best-p 0.005859375, :best-ss ["TAAACAAT" "AAATTTAC" "GTACTGTC" "TAAACGAC" "TTAACACC"]}

;; (first (sort-by #(* -1 (:best-p %)) (pmap (fn [_] (gibbs seqs 8 50)) (range 100))))
;; {:best-p 0.0087890625, :best-ss ["GTAAACAA" "CCTCGCAA" "GTACTGTC" "GTAAACGA" "CTTAACAC"]


(defn test-data [n & [bkgd]]
  (let [data-str (slurp (str "../data/data" n ".txt"))
        lines    (str/split-lines data-str)]
    (sort-by #(* -1 (:best-p %))
             (pmap (fn [_] (gibbs lines 10 50 bkgd))
                   (range 20)))))

(test-data 1)
;; below are the top results from running test-data twice
;; {:best-p 1.0, :best-ss ["ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC"]}
;; {:best-p 1.0, :best-ss ["TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC"]}
;; {:best-p 1.0, :best-ss ["AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT"]}
;;
;; From these results, we can infer a hidden motif of length 13: AATTCGAATTCCA
;; Let's make sure this still works using relative entropy as the scoring metric:

(def uniform-bases {\A 0.25 \T 0.25 \C 0.25 \G 0.25})
(test-data 1 uniform-bases)
;; Now :best-p corresponds to the entropy
;; {:best-p 20.0, :best-ss ["TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC"]}
;; {:best-p 20.0, :best-ss ["ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC"]}
;; {:best-p 20.0, :best-ss ["AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT"]}
;; Yep, we're still finding the best motif candidates!

;; Moving on to data2:
(test-data 2)

;; a sampling of some of the best results
(def test2-results [{:best-p 0.1285093120212051, :best-ss ["TCTGTCTGCT" "ACTGTCTACT" "AAAGTTTACC" "TCGGTCTACC" "TAAGTCTACT" "TCAGTCTACT" "TCTGTCTGCT" "TCTATCTACT" "TCGGTCTACT" "TCTGTCTACT"]}
{:best-p 0.11989207715519623, :best-ss ["CTGTCTATTA" "CTGTCTACTA" "CTACCAGTAA" "CTATCTAATA" "CTGCTTAATA" "CAGTCTACTA" "CTGTCAACTA" "CTATCTACTA" "CTGTCAACTA" "CTGTCTACTA"]}
{:best-p 0.0932493949309581, :best-ss ["TGTCTACTAA" "TGTCTACTAA" "TGTATACTAA" "GGTCTACCGA" "TGTCAACTAG" "CATCAGTTTT" "GGTCAACTAA" "TATCTACTAA" "TGTCAACTAA" "TGTCTACTAG"]}
{:best-p 0.08841423624506653, :best-ss ["TGTCTACTAA" "TGTCTACTAA" "TGTATACTAA" "AGCGTTACTC" "TGTCAACTAG" "AGTCTACTAC" "TGTCAACTAC" "TATCTACTAA" "TGTCAACTAA" "TGTCTACTAG"]}
{:best-p 0.07081125326122774, :best-ss ["TTGTCTACTA" "CTGTCTACTA" "TCTTCTACTA" "CTGATTACTA" "TCGGCTACTA" "CAGTCTACTA" "CTGTCTGCTA" "TAGTTTACTA" "TCTACTACTA" "CTGTCTACTA"]}
{:best-p 0.053545553218207156, :best-ss ["TCTACTAATC" "TGTACTACTC" "TCTTCTACTA" "ACAGCGTACC" "TCTACTACTC" "TCTACTACTC" "TCAACTACAA" "TCTACTAAAC" "TCTACTACTA" "TTGACACATC"]} {:best-p 0.04589618666528186, :best-ss ["ATCTGTCTGC" "ATCGGTGTAC" "ATCTCTGGGT" "ATCTGTCGAA" "ATTTCTATAA" "ATCTGTGTAC" "ATCTGTCAAC" "ATTGGGAAAC" "ATCGGTCTAC" "ATCTGTCTAC"]}
{:best-p 0.0458961842716093, :best-ss ["TTGTCTACTA" "CTGTCTACTA" "CCGTCTTCTA" "CGGTCTACCG" "AAGTCTACTA" "CAGTCTACTA" "AGCTCTTCCG" "TAGTTTACTA" "CGGTCTACTA" "CTGTCTACTA"]}
{:best-p 0.045329567211914217, :best-ss ["GTCTGCTATA" "GTCTACTACC" "GTCTTCTACT" "ATCTAATATC" "GTCTACTACT" "GTCTACTACT" "GTCTGCTATC" "ATTCATACTT" "GTCTACTACT" "GTCCTCTACA"]}
{:best-p 0.045329567211914217, :best-ss ["AATCTGTCTA" "GAGCTGTCTA" "TATATGTATA" "TATCTGTCGA" "AAGTTCCGGT" "TATCTGTGTA" "TATCTGTCAA" "GATCTATCTA" "GCTCTGTAGA" "TATCTGTCTA"]}
{:best-p 0.04371065432982424, :best-ss ["TCTGTCTGCT" "ACTGTCTACT" "ACTATCGACA" "TCTATCTAAT" "TCTATTTAGA" "ACTGACTAGA" "TCTGTCTGCT" "TCTATCTACT" "ACTTACTAAT" "TCTGTCTACT"]}])

(defn most-frequent-results [test-results]
  (take 5 (sort-by (fn [[k v]] (* -1 v))
                   (frequencies (flatten (map #(:best-ss %) test-results))))))

(most-frequent-results test2-results)

;; The 5 most commonly occuring motifs out of this sample are:
;; ["CTGTCTACTA" 6]
;; ["TGTCTACTAA" 4]
;; ["TCTGTCTGCT" 4]
;; ["CAGTCTACTA" 3]
;; ["GTCTACTACT" 3]

(test-data 2 uniform-bases)
;; Here are the top results from one round of test-data weighed by relative-entropy:

(def test2-entropy-results
  [{:best-p 15.583330452654646, :best-ss ["CTGTCTGCTA" "CTGTCTACTA" "ATATTTTTCA" "CTGTCGAATA" "GTGTCAACTA" "GTGGCTGCTA" "CTGTCTGCTA" "CTATCTACTA" "CTGTCAACTA" "CTGTCTACTA"]}
   {:best-p 14.916663766119763, :best-ss ["AATCTGTCTA" "GAACTGTCTA" "CCGGTACAAG" "CATCTATCTA" "CATGTGTCAA" "TATCTGTGTA" "TATCTGTCAA" "GATCTATCTA" "CATCTGTCAA" "TATCTGTCTA"]}
   {:best-p 14.465614508060504, :best-ss ["CTGTCTGCTA" "CTGTCTACTA" "CAAGCAACTC" "CTATCTAATA" "CAGGCTAGGC" "CTGTGTACTA" "CTGTCAACTA" "CTATCTACTA" "CTGTCAACTA" "CTGTCTACTA"]}
   {:best-p 14.251776363186782, :best-ss ["GACTACTAGA" "CTTTACTATT" "GTCTTCTACT" "GATTACTAGG" "GTCTACTACT" "GTCTACTACT" "TGTCGTTAAG" "GTTTACTACG" "GTCTACTACT" "GTCTCCTACG"]}
   {:best-p 12.326659007354099, :best-ss ["CTACTAATCT" "TTACTATTCT" "TCACAATACT" "GAACTAGACT" "TTAATAAACA" "GAAATAATCC" "CTACTAGACA" "CTACTAAACT" "TTACTAATCT" "ACAGATGTGT"]}
   {:best-p 11.757067490969789, :best-ss ["TAGGTTGGAA" "TCCTTAACGC" "CATTAAAAGC" "TCTTTCAAGC" "CACTTCAGTC" "TCCTTCAAGC" "TATGTCGAGC" "TCTTAAGAGC" "TCTTAAGGGC" "CCCTACGCGC"]}
   {:best-p 11.6710886530606, :best-ss ["ATCCAGTAAA" "GTCTACTAAG" "GTATACTAAA" "GAATCCTAAG" "TCCGCTGCTA" "ATATCTTATA" "GTATCGTAAA" "ATCTCTTAAG" "GGATACTATG" "AAACCTTATG"]}
   {:best-p 11.655648010886544, :best-ss ["ATAAATTCAT" "CCACTGGTAG" "ACAATACTAG" "ACAAGCCCGG" "AAACTGGTGG" "ACAAGCGTAG" "ACACGGGCAG" "AAACGAGTAG" "ACACTGCAAG" "AATCGCCTAG"]}
   {:best-p 11.306157659314447, :best-ss ["GCGAGCAATC" "GCGAAGTCTC" "ACGATAAGTC" "GCGAAGACTT" "GCGGTATGTC" "TAACTAGGCA" "ACCATAACTT" "GCGAACAGTC" "GCCATCTGTC" "GCGGGGTATC"]}
    {:best-p 11.256356039875568, :best-ss ["CAGCTAGATT" "ATTTCCTTAA" "ATACTAAAAC" "CTACGACAAC" "CTACTCACTT" "CTGCTAACAT" "CTGCTACCAT" "AAGCGAACAG" "CAACTAACAG" "CTACTAGCTG"]}])

(most-frequent-results test2-entropy-results)
;; Again, we see a similar set of candidate motifs:
;; ["CTGTCTACTA" 4]
;; ["CTGTCAACTA" 3]
;; ["CTGTCTGCTA" 3]
;; ["GTCTACTACT" 3]
;; ["CTATCTACTA" 2]

;; Now simply trying to maximize P is insufficient.
;; The results are skewed in favor of A/T repeats

(test-data 3)
;; {:best-p 0.13362902483849626, :best-ss ["AAAAGAAAAA" "AAAAAAAAAA" "AGAGGAAGAA" "AAAAAATAAA" "AAAAGAAAAT" "AACAGAAAAT" "AAAAATGAAA" "AAAAAAAAAA" "AGAAAAGAAA" "AAAAATTTAA" "AAAAAAAAAA" "AAAAAAAAAG" "AAAAGAGAAG" "AAAGAAAGAG" "AGAAGAAAAG" "AAAGATTATA" "AAAAATTTAA" "AAAAAAATAA" "AAAAAAAAAG" "AAAAATAGTA" "AAAAAAAATA" "AAAAAAAAAA" "AAAAAAAGAA" "AAAAGAATAA" "AAAAAAAAAA" "AAAAAAGAAG" "AAAAGAAAAA" "AAAAATATAA" "AAAGAAGGAA" "AGCAAAGAAT" "AAAAGAAAAA" "AAAAAAAAAA" "AGAAGAAAAA" "AAAAAAAAAA"]}
;; {:best-p 0.1148695910789486, :best-ss ["AAAAAAGAAA" "AAAAAAAAAA" "AAAGGAGGAA" "AAAAAATAAA" "AAAAAAAAAG" "TAAGAATAAA" "TAAAAATGAA" "AAAAAAAAAA" "TAAAGATAAA" "AAGTAATCCG" "TAAAAAAAAA" "AAAAACAAAA" "AAAAGAGAAG" "AAAGAAAGAG" "AGAAGAAAAG" "AAAGGGAGAA" "AAATGAAAAA" "AAAAAATAAA" "AAAAAAAAAG" "AAATAGTAAA" "AAAAAAAATA" "AAAAAAAAAA" "AAAAAAAGAA" "TGATAAAAAA" "AAAAAAAAAA" "AAAAAAGAAG" "AAAAAGAAAA" "AAGAAGAAAA" "AAAGAAGGAA" "AAAAACCAAA" "AAAAAGAAAA" "AAAAAAGAAA" "AAAAGAAGAA" "AAAAAAAAAA"]}
;; {:best-p 0.11364343620251935, :best-ss ["AGGAAAAAAT" "AAAAAAAAAA" "AGGAACAAAT" "AGAAAAAATA" "AAAAAAAAGA" "AAGAATAAAA" "AAAATGAAAT" "AAAAAAAAAA" "AGAAAAGAAA" "TGAAAAATTT" "AAAAAAATAA" "AAAATAATAA" "AGAAAAGAGA" "AAGAATATAA" "AGAATAAAAA" "AGAATAGTTA" "AAAATGAAAA" "AAAAAAATAA" "AAAAAAAAAA" "AAGAATAAAA" "AAAAAAAATA" "AAAAAAAAAA" "AAAATAAATA" "TAAAAAGAAT" "AAAAAAAAAA" "AAAAAGAAGA" "AAAAAGAAAA" "AAAAATATAA" "AAAATAGTAA" "AATATAAATA" "AAAAAGAAAA" "AAAAAAAAAA" "TAGAAAAAAA" "AAAAAAAAAA"]}
;; {:best-p 0.11359138357654087, :best-ss ["AAAAAAAGAA" "AAAAAAAAAG" "AAGGATAGAA" "AAATAAATAA" "CAAAAAAAAA" "AAGAATAAAA" "AAAAATGAAA" "GAAAAAAAAA" "AAGATAAAAA" "CCATGAAAAA" "AAAAAAAAAA" "CAATAAAAAA" "AAGAAAAGAG" "AAGAATATAA" "CAGAATAAAA" "AAGTAAGAAA" "CAAAATGAAA" "AAAAATAAAA" "AAAAAAAAAA" "AAGAATAAAA" "AAAACAAAAA" "CAAAAAAAAA" "AAGAAAAAAA" "AAGAATAAAA" "AAAACAAAAA" "GAAAAAAGAA" "GAAAAAGAAA" "AAAAATATAA" "AAAATTGAAA" "CGATTAAAAA" "AAAAAAGAAA" "AAAAAAAAAA" "AAGGAAAAAA" "AAAAAAAAAA"]}

;; No good!

(def low-gc-dist {\A 0.315 \T 0.315 \C 0.185 \G 0.185})
(test-data 3 low-gc-dist)
;; Still doesn't work well:
;; {:best-p 11.52886007980086, :best-ss ["AAAAAAAGAA" "AAAAAAAGAA" "ATGAGTGGAA" "AAAAATAGAA" "ATAAAAATAA" "AAGAATAAAA" "AAAAATGAAA" "AAAAAAAAAA" "AAGAAAAGAA" "AAAAATTTAA" "AAAAATAAAA" "AAAAAAAAAG" "AAGAAAAGAG" "AAGAATATAA" "AAGAAAAGAG" "AAGTAAGAAA" "AAAAATTTAA" "AAAAAATAAA" "AAAAAAAAAA" "AAGAATAAAA" "AAAAAAATAG" "AAAAAAAAAA" "AAAAAAAGAA" "AAAAGAATAA" "AAAAAAAAAA" "AAAAAAGAAG" "AAGAAAAAAG" "AAAAATATAA" "AAAAATTGAA" "ATAAATATAA" "AAAAAAGAAA" "AAAAAAAAAA" "AAAAGAAGAA" "AAAAAAAAAA"]}
;; {:best-p 11.371393724675492, :best-ss ["AAAAAAGAAA" "AAAAAAAAGA" "AAAGGATAGA" "AAAAAATAAA" "AAAAAAAAGA" "AGATAAGAAA" "AAAAAATTAA" "AGAAAAAAGA" "AGAAAAGAAA" "CCATGAAAAA" "AAAAAATAAA" "AGAAAAATAA" "AAAGAAAAGA" "AAAAGAAAGA" "AGAAAAGAGA" "AATGGAAAAA" "AAATGAAAAA" "AAAAAATAAA" "AAAAAAAAGA" "AGTAAAAAGA" "AAAACAAAAA" "AAAAGAATAA" "AAAAAAATAA" "AAAAGAATAA" "AAAAAAAAAA" "AGAAGAGAAA" "AAAGAAAAAA" "AGAAGAAAAA" "AAAGAAGGAA" "AAAGCAAAGA" "AAAAGAAAAA" "AAAAAAAAAA" "AAAGAAGAAA" "AAAAAAAAAA"]}
;; {:best-p 11.327796258009817, :best-ss ["AAAAAGAAAA" "AAAAAAAAAA" "AAAGTGAGAA" "GAAAAAATAA" "AAAAAGATAA" "AAGAATAAAA" "AAAAAAATTA" "AAAAAAAAAA" "AAGATAAAAA" "GAAAAATTTA" "AAAAAAATAA" "AATAAAAAAA" "GAAAAGAGAA" "AAGAATATAA" "GAAAAGAGAA" "AAAGATTATA" "AAAATGAAAA" "AAAATAAAAA" "AAAAAAAAAA" "AAGAATAAAA" "AAAAAAAATA" "AAAAAAAAAA" "AAAAAAAGAA" "AAAAAGAATA" "AAAAAAAAAA" "GAAAAAAGAA" "AAAAAGAAAA" "AAGAAGAAAA" "AAAAATAGTA" "AATATAAATA" "AAAAAGAAAA" "AAAAAGAAAA" "GAAAAAAGAA" "AAAAAAAAAA"]}
;; {:best-p 11.238083345541527, :best-ss ["AAAAGTAAAA" "AAAAAAAAAA" "TAGCGAAAAA" "TAGAGAAAAA" "ACAAAAAAAA" "AAGAATAAAA" "AAGAGAACAA" "AAAAAAAAAA" "AAGATAAAAA" "AAAAATTTAA" "AAAAAAAAAA" "TAAAAAAAAA" "TCAAAGAAAA" "AAGAATATAA" "CAGAATAAAA" "CAGAGAAAAA" "TATAATAAAA" "AAAATAAAAA" "AAAAAAAAAA" "AAGAATAAAA" "AAAAAAAATA" "AAAAAAAAAA" "AAGAAAAAAA" "AAGAATAAAA" "TAAAAAAAAA" "AAGAGAAATA" "AAAAGAAAAA" "AAGAAGAAAA" "AAACGTTTAA" "AAGAATACAA" "AAAAGAAAAA" "AAAAAAAAAA" "TAGAAAAAAA" "AAAAAAAAAA"]}


(test-data 4 low-gc-dist)
;; {:best-p 10.131471123522752, :best-ss ["TTTTTTTTTC" "TTTTTGTTTT" "TTTTTGCTTC" "TTTTTTTTTC" "TTTTTTTTTT" "TATTCTTTTT" "CTTTTTTTTT" "TTTTCTTTTT" "TACTCGTTTT" "TATTTTTTTT" "TTTTTTTTTT" "TTTTTACTTT" "CTTTTTTTTC" "TTCTTTTTTT" "TTTTTTGTTG" "CTTTTTCTAT" "TTTTTTTTAT" "CATTTATTTC" "TTTTTTTTTT" "TTATTTTTTT" "TTTTTTTTTT" "TACTTTCTTT" "CTTTTATTGT" "ATATCTTTTC" "TTTTCTTTTT" "TTTTTTCTTC" "TTCTTTCTAT" "TTTTTTCTTT" "TTTTCTGTTT" "TTTTTTTTTT" "CTATTTTTTC" "TTTTTTCTTC" "CTTTTTTTTC" "TTTTTATTTT" "TTTTTTCTTC" "TTTTCATTAT" "TATATTCTTG" "TACTTTCTTT" "TATTTTTTTC" "TTTTTTCTTC" "TTCACTCTTT" "TTTTCGGTTT" "TTTACTTTTG" "TATTTGCTTT" "TTTTTTTTTT" "CTTTTGTTTT" "CTATTGTTTC" "TTTACTTTTT" "TTAACTCTTG" "TTTACTTTTT" "TATTTGCTTT" "TTCTTATTTC" "TTTTCTTTTC" "TTTTTTCTAT" "TTTTTTTTTT" "TTTTTTTTTT"]}
;; {:best-p 10.040892286282462, :best-ss ["ATATATATAT" "ATATATCAAA" "AAAAGACAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAATGAAAAA" "AACAAAAAAA" "ATATATACAT" "AAAAGAAAAA" "AAATAAACAA" "AAAAAAAAAA" "AAAAGAAAAT" "AAAAAAAAAA" "AAAAGAAAAA" "AAGAAATAAT" "AAGAATAAAA" "AAATTAAAAT" "AAAAAAAAAA" "AACAAAAAAA" "AAAAAAAAAA" "AAATGAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAG" "ATAAATAAAT" "AAAAAAAAAA" "AAAAAAAATA" "AACTGTTTAA" "AAAAAAAAAT" "ATAAAAAAAA" "AAAAGATATA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAATGAAATA" "ATAAATAAAT" "AAAAAAAATA" "AAATAAAAAG" "AACAAATAAA" "AAGAAAAAAT" "ATATATAAAA" "ATATATAAAG" "AAATATAAAT" "AAATTTAAAA" "ATATAAAAAA" "AAAAAAAAAG" "AAATAAATAA" "AAAAATTAAA" "AAAAAATAAA" "ATAAGCAATA" "AACAAAAAAA" "AAAAAAAAAG" "AACAATAAAA" "AAATACAATT" "AAAAATAAAA" "AAATACAAAA"]}
;; {:best-p 10.021928335342775, :best-ss ["TAAATTAAAA" "AACATATATC" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AATAAAAACA" "AAAAAAAACC" "TATATATATA" "AAAAGAAAAA" "AACAAATAAA" "AAAAAAAAAA" "AATAAAAAAA" "AAAAAAAAAA" "AAAAGAAAAA" "AATAGATAAC" "TATATATAAA" "TAAATTAAAA" "AAAAAAAAAA" "AAAAAATAAC" "AAAAAAAAAA" "TATAATAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AAAAAAAAAA" "AATAAATACA" "AAAAAAAAAA" "AAAAAAAATA" "CAAAAATAGA" "TAAAAAAAAA" "AAAAGAAAAA" "AAAAGATATA" "AAAAAAAAAA" "AAAAAAAAAA" "CATATTAATA" "ATAAAAAATA" "TAAATAAATA" "AAAAAAAATA" "ATTAAAAAAA" "AAAAAAAAGA" "AAAAAATATC" "TATATAAAAC" "CATATATAAA" "AATATAAATA" "AAAATTTATA" "AAAAATTAAA" "AAAATAAATA" "AAAATAAATA" "AAAATTAAAA" "AAAAATAAAA" "TATAATAAGC" "AAAATAAACA" "AAAAAAAAAA" "AATATTAAAA" "TACAATTACA" "AAAAAAAACA" "AATAAAAACA"]}

;; This faux-relative entropy approach is obviously not working well.
;; Going to try weighted discarding based on relative entropy now.