cat ~/athena/ukbiobank/hesin/hesin.txt | fgrep -w -f search_eid.txt > small_hesin
cat ~/athena/ukbiobank/hesin/hesin_diag.txt | fgrep -w -f search_eid.txt > small_hesin_diag
cat ~/athena/ukbiobank/hesin/hesin_oper.txt | fgrep -w -f search_eid.txt > small_hesin_oper
cat ~/athena/ukbiobank/hesin/death.txt | fgrep -w -f search_eid.txt > small_death
