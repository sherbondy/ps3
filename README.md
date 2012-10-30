I decided to implement the gibbs sampler in clojure. See: `code/gibbs/src/gibbs/core.clj`

You'll probably want to run the test-data function:
```
(test-data 1 10 2)
```

Do `lein repl` and enter the above command.
If you don't have Leiningen, install it by following the directions here:
[http://leiningen.org/](http://leiningen.org/)

The above example seeds 2 gibbs-samplers and runs them in parallel to find the globally most likely 10-mer motifs for data1.txt.


