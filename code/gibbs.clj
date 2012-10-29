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
(def uniform-bases {\A 0.25  \T 0.25  \C 0.25  \G 0.25})  ; data 1 and 2
(def low-gc-dist   {\A 0.315 \T 0.315 \C 0.185 \G 0.185}) ; data 3 and 4

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
;; (reduce + (vals (pfreq-dist "ATGCC")))

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
    (min (count (filter #(< % rand-f) prob-cdf))
         (dec (count prob-dist)))))

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
           ;; arbitrarily set convergence cutoff to be 20 rounds
           ;; without improvement
             (pmap (fn [_] (gibbs seqs len 20 bkgd))
                   (range reps))))

;; (first (n-gibbs seqs 8 100))
;; {:best-p 0.0098876953125, :best-ss ["GTAAACAA" "CCTCGCAA" "GTCAAGCG" "GTAAACGA" "CTTAACAC"]}

(defn test-data [data-i len & [bkgd]]
  (let [data-str (slurp (str "../data/data" data-i ".txt"))
        lines    (str/split-lines data-str)]
    (n-gibbs lines len 10 bkgd)))

;; (test-data 1 10)

;; below are the top results from running test-data twice
;; {:best-p 0.4666044141799745, :best-ss ["AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT"]}
;; {:best-p 0.3639411606231313, :best-ss ["AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "AATTCGAATT" "GATTCGAATT" "AATTCGAATT" "CATTCGAATT"]}
;; {:best-p 0.1863975964873591, :best-ss ["CGAATTCGAA" "CGAATTCGAA" "GGAATTCGAA" "CGAATTCGAA" "GCAATTCGAA" "CGAATTCGAA" "CGAATTCAAA" "GCAATTCGAA" "CCAATTCGAA" "CGAATTCGAA"]}

;; (most-frequent-results (test-data 1 10))
;;

;; Moving on to data2:
;; (test-data 2 10)

;; (most-frequent-results (test-data 2 10))

;; a sampling of some of the best results
(def test2-results
  [{:best-p 0.11931785733171743, :best-ss ["CTGTCTGCTA" "CTGTCTACTA" "CCGTCTTCTA" "CTATCTAATA" "GTGTCAACTA" "CTGTGTACTA" "CTGTCTGCTA" "CTATCTACTA" "CTGTCAACTA" "CTGTCTACTA"]} {:best-p 0.08156155203070686, :best-ss ["CTGTCTATTA" "CTGTCTACTA" "CCGTCTTCTA" "CTGTCGAATA" "GTGTCAACTA" "CTGTGTACTA" "CTGTCTGCTA" "CTATCTACTA" "CGGTGTAGGA" "CTGTCTACTA"]} {:best-p 0.06125861119679491, :best-ss ["TGTCTACTAA" "TGTCTACTAC" "TGTATACTAA" "AGACTACGAA" "AGTCTACTAC" "AGTCTACTAC" "TGTCAACTAC" "TATCTACTAA" "TGTCAACTAA" "TGTCTCCTAC"]} {:best-p 0.0181841514780996, :best-ss ["ACTTCCATAC" "ATTTCCTTAA" "ATGATCGCAA" "TCTTCCGGAA" "ATGTCCACAA" "AGGTCCTTCA" "ATTTTCGCAA" "ATTTCCGTCA" "ATTTCCAGAA" "ATTTCCGTAA"]} {:best-p 0.004044962046614947, :best-ss ["ACCATAGACT" "AGCTATGACT" "ACAATCGACG" "AGAACAGACT" "AGTCTAAACT" "AGTTTTGACT" "AGAAAAGACT" "AGTCGCGACT" "AGAATTGACT" "CAATGAAACT"]} {:best-p 0.003770681847330107, :best-ss ["CCATACGGCG" "ATGAACTGCG" "CCGAAGGGGA" "CTTAACAGGA" "CTTTACTGCA" "CTGTACACGG" "CTGTACGGCG" "CTTTCCGGCA" "ATGAAGAGCA" "CCGAACGCGG"]} {:best-p 0.0036435117343088306, :best-ss ["CAGCTAGATT" "CCGAAAGCGT" "CAGTAATCTC" "CGCTAAGCTC" "CGGAAGTCTT" "CGCATATCTT" "TCCTTGGCTC" "ACGTAATCTC" "CGGCACGCTG" "CAGCACGCGC"]} {:best-p 0.0010550766281002353, :best-ss ["GCTTCTAGCC" "GCGGATAACA" "CCAGCTTTCA" "ACGGCTTGCC" "GCTGCTAAGC" "GCTGCTAACA" "CCGGCTATGA" "ACGCCTTGGA" "GCGGCCAGGC" "CCGGCCTTGA"]} {:best-p 9.445888599080243E-4, :best-ss ["GAGCAATCAC" "GCGCAACGGA" "TCGCAACGGC" "TCTCAACGGC" "GCGCCATGGC" "TAGGCAGCGA" "TCCCAGCCGC" "GCCGTGAGGG" "GCCGCAACAC" "TACGCGCCGA"]} {:best-p 3.1102138955200765E-4, :best-ss ["CATCCCAGCT" "TGAGCCACGT" "CATCCCCGTT" "TAAGCTCCCC" "CGAGCTACTT" "CGAGTCAGTG" "TATGCGGGTT" "TCATTTCCGT" "CATGTCACGC" "TAACTCCCGT"]}])

;; (most-frequent-results test2-results)

;; The 5 most commonly occuring motifs out of this sample are:
;; ["CTGTCTACTA" 4]
;; ["CTGTCTGCTA" 3]
;; ["CTATCTACTA" 2]
;; ["GTGTCAACTA" 2]
;; ["AGTCTACTAC" 2]

;; another round of results:
;; (["ACTCTACCAG" 2] ["CTGTCTGCTA" 2] ["ATCTGTCTGC" 2] ["CTGTCTACTA" 2] ["CTACTAGAAG" 1])


;; For 3 and 4 simply trying to maximize P is insufficient.
;; The results are skewed in favor of A/T repeats

;; (test-data 3 10 low-gc-dist)

;; Still doesn't work very well:
;; {:best-p 0.009131122250906587, :best-ss ["AAAAGTAAAA" "AAAAAAAAGA" "AAAAGGAGGA" "AAAAAATAAA" "AAAAAAAAGA" "AATAAAAAAC" "AAAAAAATTA" "AAAGAAAAAA" "AAGAAAAGAA" "AAAAATTTAA" "AAAAAAAAAA" "AATAAAAAAA" "AAAGAAAAGA" "AAAAGAAAGA" "TAGAAGAAAA" "CAGAGAAAAA" "AAAATGAAAA" "AAAAATAAAA" "AATAAAATAA" "AATAGTAAAA" "AAAAAAAATA" "AAAAAAAAAA" "AAGAAAAAAA" "AAAAAGAATA" "AAAAAAAAAA" "AAAAAGAAGA" "AAAAAGAAAA" "GAAAAAAATA" "AAAATAAGTA" "AAGCAAAGAA" "AAAAAGAAAA" "AAAAAAAAGA" "TAGAAAAAAA" "AAAAAAAAAA"]}
;; {:best-p 0.006884564319836377, :best-ss ["AAAAAAGAAA" "AAAAAAAAAA" "AGAAAATACC" "AAAAAATAAA" "AAAAAAAAGA" "AGAATAAAAA" "AACAACAAAT" "AAAAAAAAAA" "AGAAAAGAAA" "AAAATTTAAC" "AAAAAAAAAA" "AAAAACAAAA" "AGAAAAGAGA" "ATATAAGAAT" "AATAACAAAA" "AAAAATGCAA" "AAAATCTCAA" "AAAAAATAAA" "AAAAAAGAGA" "AGAATAAAAG" "AAAACAAAAA" "AAAAAAAAAA" "AAAAAATAAC" "AAGAATAAAA" "AAAAAAAAAA" "AAAAAGAAGA" "AAAAGAAAAA" "AAAATATAAT" "AAAATTGAAA" "AAATATAAAT" "ATAATATAAA" "AAAAAAAAAA" "AAAAAAAGAA" "AAAAAAAGAA"]}
;; {:best-p 0.004072147475045898, :best-ss ["TTTTTATTTA" "TTTTTCTTTT" "TTTTTCATTT" "TTTTTTTTTT" "GATTTTTTTT" "TTTCTATTAT" "TCTCTTCCAT" "TTATTCCCTT" "TTGTATATAT" "TTGCTCAAAA" "TCTTTATTAT" "ACTTTTTTTT" "TTTTGCCTAT" "TTTGTCTTTA" "TTTTTATTTT" "TTTTTTTTTT" "TTTTTTTTTT" "TTTTTCTTAT" "TTATTATTTT" "TTTTTTTTTT" "TTTCTTCTTT" "TTTTACTTAA" "TTTTTAATAT" "TTTTTTTCTT" "TTTTTTTCTT" "TTTCTAGTAT" "TTTTTCCTTT" "TTTTTTTCTT" "TTGTTTCTTA" "TTTTACTTTT" "TTTTTTTTTT" "TTTTGTTTTT" "GATCTCTTTT" "TTTCTTTTTT"]}
;; {:best-p 0.002952289977243015, :best-ss ["AAAAAAAACA" "AGAAAAAGAA" "AGAAAATACC" "GAAAAAATAA" "AAAAAAAAGA" "ATAAAAAACA" "AAAAAAATTA" "AGAAAAAAGA" "AAGAAAAGAA" "TGAAAAATTT" "AAAAAAAACA" "AAAAAAAACA" "AGAAATAAGG" "AGAAATATTA" "ATAAAAAACA" "AGAAAAACCA" "TTAAATAATA" "TAAAAAAATA" "AAAAAAAAGA" "AAAAATAGTA" "AAAAAAAATA" "AAAACAAGCA" "AGAAAAAGGA" "TTAAAAAGAA" "GGGAAAAAAA" "TGAAAAAAGA" "TGGAAAAGCC" "AGAAAAAATA" "AAAAATAGTA" "TAAAAAACCA" "AAAAAAAACA" "AAAACAAACC" "AAAAAAAGAA" "AAAAAAAAAA"]}
;; {:best-p 0.002947230754676865, :best-ss ["TAATAAAAAT" "GGATAAAAAA" "GATAATACAA" "CAGAAAAAAT" "GGAAAAAAAA" "GAATAAAAAA" "GAGAAAAAAA" "CAGAAAAAAA" "GAAATTAAAA" "CACAGAATAA" "AAAAAAAAAA" "GATATTAAAA" "TAGAAATAAG" "AAAAGAAAGA" "GAATAAAAAA" "CAGAGAAAAA" "CAGTTTAAAA" "CTTAAAAAAA" "GAAAAAAAAA" "TTAAAAATAG" "GAAAAAAAAT" "TAGAATTAAA" "GAAAAAAAGA" "AAAAGAATAA" "AAAAAAAAAA" "GACAATTCAA" "AAAAGAAAAA" "TAGAAAAAAA" "AAAATAGTAA" "GAATGTATAA" "CAAAAAAAAA" "TCTAAAAAAA" "AGAAAAAAAG" "GAAAAAAAAA"]}
;; {:best-p 4.3920720567776453E-4, :best-ss ["TTTTTCATCC" "TTTTTCGCTG" "GTTTATATTT" "CTTTTTTTTT" "TTTTGCAAGT" "TGCATAATTG" "ATTTTTTTTG" "TTATTTATGT" "TTGTATATAT" "TTCTTATGTT" "CTCTTTGGAA" "CTTTTTATTA" "TTCAAATTTC" "TTGTTTGTCT" "TTTTATTTTA" "TTTTATTTTT" "TTTTTTTTTT" "GCTTATAAAT" "TTTTCTATTA" "TTTTTTTTTT" "CTGTACATAT" "TTGAATGGTA" "TTTTAATATA" "TTGTATATAT" "CTTTAAATAT" "CTTTACATAT" "CTTTCTGCAC" "GTGTGCGTGA" "TTTTTAGTTT" "CTTAAATATA" "TTTTTTTATT" "TTTTGTTTTT" "ATTATTTATT" "TTTACTAAAT"]}

;; most-frequent-results on this data yields:
;; (["AAAAAAAAAA" 14] ["TTTTTTTTTT" 9] ["AAAAAAAAGA" 6] ["AAAAAAAACA" 4] ["TTGTATATAT" 3])

;; Out of these, only "TTGTATATAT" is likely to be a genuine motif.

;; (most-frequent-results (test-data 4 10 low-gc-dist))

;; ["AAAAAAAAAA" 34] ["TTTTTTTTTT" 21] ["AAAAAAAAAT" 5] ["AACAAAAAAA" 5] ["AAAAATTAAA" 4]

;; ["AAAAAAAAAA" 42] ["AACAAAAAAA" 8] ["TTTTTTTTTT" 5] ["AAAAAAAAAT" 5] ["ATAAATAAAT" 5]