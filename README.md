# Parallel-and-Distributed-Systems  
### Κώδικας 1ης Εργασίας  
Παυλίδης Χρήστος ΑΕΜ: 9480 pavlidic@ece.auth.gr

---

Τα mmio, coo2csc, timediff και Μ312 περιέχουν κάποιες συναρτήσεις που τα υπόλοιπα αρχεία χρησιμοποιούν για:
- ανάγνωση πινάκων
- μετατροπή από coo σε csc μορφή
- χρονομέτρηση και
- πολλαπλασιασμό πινάκων 

---

Οι υλοποιήσεις με cilk γίνονται compile με clang καθώς με gcc γίνονται πολύ πιο αργές.  

Παράδειγμα:
`{clang bin} -o {output name} -O3 {.c file} mmio.c coo2csc.c -fcilkplus`
