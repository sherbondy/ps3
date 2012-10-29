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
(def uniform-bases {\A 0.25 \T 0.25 \C 0.25 \G 0.25})

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

(defn pfreq-dist
  "Freq distribution with pseudocounts"
  [v]
  (let [o-dist (freq-dist v)
        pseudo 0.01
        fcount (float (count v))]
    (into {}
          (for [base bases]
            (let [o-freq (get o-dist base 0.0)]
              {base (/ (+ o-freq pseudo)
                       (+ 1 (* pseudo fcount)))})))))

;; (pfreq-dist "TCCGGCCC")
;; ((freq-dist "ATGCC") (nth bases 1))
;; (reduce + (vals (freq-dist "ATGCC")))
;; (freq-dist (map #(nth % 0) (subseqs seqs 8)))

(defn make-profile [subseqs size bases]
  ;; first index = letter
  (let [arr (make-array Float/TYPE 4 size)
        bcount (count bases)]
    (doseq [i (range size)]
      (let [freqs (pfreq-dist (map #(nth % i) subseqs))]
        (doseq [j (range bcount)]
          (aset-float arr j i
                      (get freqs (nth bases j)))
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

(defn lmer-prob [lmer profile bkgd]
  (reduce * (map-indexed (fn [i base]
                           (get-in profile [(bi base) i]))
                         lmer)))

(defn log2 [n]
  (/ (Math/log n) (Math/log 2)))

(defn rel-entropy [len profile bkgd]
  (reduce + (for [j (range len)]
              (reduce + (for [r bases]
                          (let [p-rj (get-in profile [(bi r) j])]
                            (* p-rj (log2 (/ p-rj (bkgd r))))))))))
                            
;; from http://rosettacode.org/wiki/Probabilistic_choice#Clojure
(defn to-cdf [pdf]
  (reduce
    (fn [acc n] (conj acc (+ (or (last acc) 0) n)))
    []
    pdf))

; (to-cdf [1 2.0 1 1])

(defn rand-cdf-idx
  "takes a distribution and returns a random
   index weighed by the cdf of the distribution.
   The dist does NOT have to sum to one."
  [prob-dist]
  (let [prob-cdf   (to-cdf prob-dist)
        rand-f     (rand (last prob-cdf))]
    (count (filter #(< % rand-f) prob-cdf))))

(defn rank-lmers [lmers profile bkgd]
  (let [probs      (map #(lmer-prob % profile bkgd) lmers)
        nonzeros   (filter #(not= % 0.0) probs)
        max-prob   (apply max probs)
        min-prob   (if (not (zero? (count nonzeros)))
                     (apply min nonzeros))
        prob-dist  (if min-prob
                     (map #(/ % min-prob) probs)
                     (repeat (count lmers) 1))
        new-start  (rand-cdf-idx prob-dist)]
    [new-start max-prob]))

;; this should yield a uniform dist over start positions
(defn test-ranker []
  (rank-lmers ["AAA" "ATG" "CAG"]
              (to-array-2d [[0.5 0.5 0.5]
                            [0.0 0.5 0.0]
                            [0.0 0.0 0.5]
                            [0.5 0.0 0.0]]) uniform-bases))

;; it works well
(frequencies
 (for [i (range 10000)]
   (test-ranker)))

(defn rel-profile-entropies [len ss bkgd]
  (for [i (range (count ss))]
    (let [new-ss      (remove-nth ss i)
          profile     (make-profile new-ss len bases)]
      (rel-entropy len profile bkgd))))

(defn rand-entropy-idx
  "return a random index from the distribution of
   possible removal indices weighed by relative entropy
   of the resulting profile"
  [len ss bkgd]
  (let [rel-entropy-dist (rel-profile-entropies len ss bkgd)]
    (rand-cdf-idx rel-entropy-dist)))

(defn gibbs-iter [seqs len ss bkgd]
  (let [removed-idx (rand-entropy-idx len ss bkgd)
        removed-seq (nth seqs removed-idx)
        new-ss      (remove-nth ss removed-idx)
        profile     (make-profile new-ss len bases)
        [new-start max-prob] (rank-lmers (l-mers removed-seq len)
                                         profile
                                         bkgd)
        added-ss    (l-mer removed-seq new-start len)]
    ;;(pprint new-ss)
    ;;(pprint profile)
    ;;(pprint max-prob)

    ;; must return new-ss in order!
    ;; ss must be a vector for this to work
    [max-prob
     added-ss
     (reduce conj (conj (subvec ss 0 removed-idx)
                        added-ss)
             (subvec ss (inc removed-idx)))]))

;; Idea: stop iterating once we've gone cutoff rounds
;; without seeing an increase in max-prob
;; If no background is provided, a uniform base distribution is assumed

(defn gibbs [seqs len cutoff & [bkgd]]
  ;; (println "Beginning gibbs")
  ;; rounds = the number of rounds without improvement
  (let [bkgd (if bkgd bkgd uniform-bases)]
    (loop [ss (subseqs seqs len) best-prob 0.0 rounds 0 best-ss []]
      (let [[max-prob added-ss new-ss] (gibbs-iter seqs len ss bkgd)
            new-prob                   (max max-prob best-prob)
            best-ss                    (if (= new-prob max-prob) new-ss best-ss)
            rounds                     (if (> new-prob best-prob)
                                         0
                                         (inc rounds))]
        
        (if (< rounds cutoff)
          (recur new-ss new-prob rounds best-ss)
          {:best-p new-prob :best-ss best-ss})
        ))))

;; (gibbs seqs 8 50)
;; {:best-p 0.005859375, :best-ss ["TAAACAAT" "AAATTTAC" "GTACTGTC" "TAAACGAC" "TTAACACC"]}

(defn most-frequent-results [test-results & [n]]
  (let [n (if n n 5)]
    (take n (sort-by (fn [[k v]] (* -1 v))
                     (frequencies (flatten (map #(:best-ss %) test-results)))))))


;; Repeat gibbs with reps seeds: returns the sorted results
(defn n-gibbs [seqs len reps & [bkgd]]
  (sort-by #(* -1 (:best-p %))
             (pmap (fn [_] (gibbs seqs len 50 bkgd))
                   (range reps))))

;; (first (n-gibbs seqs 8 100))
;; {:best-p 0.0098876953125, :best-ss ["GTAAACAA" "CCTCGCAA" "GTCAAGCG" "GTAAACGA" "CTTAACAC"]}

(defn test-data [data-i len & [bkgd]]
  (let [data-str (slurp (str "../data/data" data-i ".txt"))
        lines    (str/split-lines data-str)]
    (n-gibbs lines len 20 bkgd)))

(test-data 1 10)
;; below are the top results from running test-data twice
;; {:best-p 1.0, :best-ss ["ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC" "ATTCGAATTC"]}
;; {:best-p 1.0, :best-ss ["TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC" "TTCGAATTCC"]}
;; {:best-p 1.0, :best-ss ["AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT"]}
;;
;; From these results, we can infer a hidden motif of length 13: AATTCGAATTCCA

;; Moving on to data2:
(test-data 2 10)

(most-frequent-results (test-data 2 10))

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

(most-frequent-results test2-results)

;; The 5 most commonly occuring motifs out of this sample are:
;; ["CTGTCTACTA" 6]
;; ["TGTCTACTAA" 4]
;; ["TCTGTCTGCT" 4]
;; ["CAGTCTACTA" 3]
;; ["GTCTACTACT" 3]

;; another round of results:
;; (["CTGTCTACTA" 4] ["CTGTCAACTA" 3] ["GTCTACTACT" 3] ["AGAGCTCATC" 2] ["CAGTCTACTA" 2])


;; For 3 and 4 simply trying to maximize P is insufficient.
;; The results are skewed in favor of A/T repeats

(def low-gc-dist {\A 0.315 \T 0.315 \C 0.185 \G 0.185})
(test-data 3 10 low-gc-dist)
;; Still doesn't work very well:
;; {:best-p 0.1289751937490981, :best-ss ["AAAAGAAAAA" "AAAAAAAAAA" "AAAGGATAGA" "AAAAAATAAA" "CAAAAAAAAA" "ATAATAACAA" "AAAAAATTAA" "AAAAAAAAAA" "AAAGATAAAA" "AAAAATTTAA" "AAAAAAAAAA" "ATAAAAAAAA" "AAAGAAAAGA" "AAAAGAAAGA" "ATAGAAGAAA" "AAAAATGCAA" "CAAAATGAAA" "AAAAAAATAA" "AAAAAAAAAA" "CAAGAATAAA" "AAAAAATAGA" "AAAAGAATAA" "CAAAAAGAAA" "AAAGAATAAA" "AAAAAAAAAA" "AAAGAAGAGA" "AAAGAAAAAA" "AAAAATATAA" "AAAATAGTAA" "ATAAATATAA" "AAAAAAGAAA" "AAAAAAGAAA" "ATAGAAAAAA" "AAAAAAAAAA"]}
;; {:best-p 0.11381719298139799, :best-ss ["AAAAAAGAAA" "AAAAAAAAAA" "GATAATACAA" "AAAAAATAAA" "AATAAAAAAA" "AATAAAAAAC" "AAATAATCAA" "AAAGAAAAAA" "AAAGATAAAA" "AAAATTTAAC" "AAAAAAAAAA" "AATAAAAAAA" "AAAGAAAAGA" "TATATATAAA" "ACTTAATAAC" "AAAAATGCAA" "TATAATAAAA" "AAAAAATAAA" "AAAAAAAAAA" "AAATAAGAAC" "AAAAAAAATA" "AAAAAAAAAA" "AAAAAAGAAC" "AAAGAATAAA" "AAAAAAAAAA" "AAAGAAGAGA" "AAAGAAAAAA" "TAAAAATATA" "AAAATTGAAA" "AATATAAATA" "AAAAAAGAAA" "AAAAAAAAAA" "AAAGAAGAAA" "AAAAAAAAAA"]}
;; {:best-p 0.11301851134502555, :best-ss ["AAAAAAGAAA" "AAAAAAGAAA" "CAAAGTGAGA" "AAAAAATAAA" "AAATAAAAAT" "ATAAGAATAA" "AGAAAAAAAT" "AAAAAAAAAA" "AGAAAAGAAA" "AAAAATTTAA" "AAAAAAACAA" "AAAAAAACAA" "AGAAAAGAGA" "AAAAGAAAGA" "AAAAGTAAGA" "AAAAATGCAA" "AAATGAAAAA" "AAAAAAATAA" "AAAAAAAAAA" "AATAGTAAAA" "AAAAAATAGA" "AAAAAAAAAA" "AAAAAAATAA" "AAAAGTTAAA" "AAAAAAAAAT" "AGAAGAGAAA" "AAAAGAAAAA" "AAAAATATAA" "AATAATACAT" "ATAAATATAA" "AAAAAAGAAA" "AAAAAAAAAA" "AGAAGAAAAA" "AAAAAAAAAA"]}
;; {:best-p 0.10731869081755462, :best-ss ["ATACAAAAAA" "AAAAAAAAAA" "AGGAAGAAAT" "AAAAAATAAA" "AAAAAGAAAA" "AATAAAAAAC" "AAAAAAATTA" "AAAAAGAAAA" "AAGATAAAAA" "ATGAAAAATT" "CAAAAAAAAT" "CTAAAAAAAA" "AAAAGATATC" "AAAAGAAAGA" "AAAAGAGAAC" "AGTAAGAAAC" "AAAAAGAATC" "AAAAATAAAA" "AAAAAAAAAA" "AAGAATAAAA" "AAACAAAAAT" "AAAAAAAAAA" "AAAAAGAAAC" "AAGAATAAAA" "AAAAAAAAAA" "AAAAAGAAGA" "AAAAAGAAAA" "AAAAAATATC" "AAGAAGGAAA" "AAAAACCAAA" "ACACAAAAAA" "AAAAAAAAAA" "AAAAAGAAGA" "AAAAAAAAAA"]}

;; most-frequent-results on this data yields:
;; ["AAAAAAAAAA" 78] ["TTTTTTTTTT" 19] ["AAAAAGAAAA" 19] ["AAAAAATAAA" 13] ["AAAAAAGAAA" 12] ["AAAAAAATAA" 9] ["AAAAAAAAAC" 7] ["AAGAATAAAA" 7] ["GAATAAAAAA" 7] ["AAAAATTTAA" 7]

(test-data 4 10 low-gc-dist)
;; {:best-p 10.131471123522752, :best-ss ["TTTTTTTTTC" "TTTTTGTTTT" "TTTTTGCTTC" "TTTTTTTTTC" "TTTTTTTTTT" "TATTCTTTTT" "CTTTTTTTTT" "TTTTCTTTTT" "TACTCGTTTT" "TATTTTTTTT" "TTTTTTTTTT" "TTTTTACTTT" "CTTTTTTTTC" "TTCTTTTTTT" "TTTTTTGTTG" "CTTTTTCTAT" "TTTTTTTTAT" "CATTTATTTC" "TTTTTTTTTT" "TTATTTTTTT" "TTTTTTTTTT" "TACTTTCTTT" "CTTTTATTGT" "ATATCTTTTC" "TTTTCTTTTT" "TTTTTTCTTC" "TTCTTTCTAT" "TTTTTTCTTT" "TTTTCTGTTT" "TTTTTTTTTT" "CTATTTTTTC" "TTTTTTCTTC" "CTTTTTTTTC" "TTTTTATTTT" "TTTTTTCTTC" "TTTTCATTAT" "TATATTCTTG" "TACTTTCTTT" "TATTTTTTTC" "TTTTTTCTTC" "TTCACTCTTT" "TTTTCGGTTT" "TTTACTTTTG" "TATTTGCTTT" "TTTTTTTTTT" "CTTTTGTTTT" "CTATTGTTTC" "TTTACTTTTT" "TTAACTCTTG" "TTTACTTTTT" "TATTTGCTTT" "TTCTTATTTC" "TTTTCTTTTC" "TTTTTTCTAT" "TTTTTTTTTT" "TTTTTTTTTT"]}


;; Going to try weighted discarding based on relative entropy now.