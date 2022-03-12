# unattached_phd_work

This repository contains all of the projects that I have not deposited into any other repository. There are many, many scripts that are admitedly not well documented but do aim to define the extent of my work and breadth of knowledge.  The specific projects include:

**SPORE** - Creation of ancestry-aware polygenic risk scores that predict the risk of prostate cancer.  The methods employed, are hoped/hypothesized to work well on all individuals regardless of population, include IMPACT, XPASS, PolyFun and PRS-CSx.

**covid_whole_exome** - In the beginning of the pandemic surgery residents collected samples from COVID patients (and clinicians treating those patients) that were then whole exome sequenced. I attempted to determine if any variants associate to COVID or its outcomes, but the sample size was far too small.

**exome_score** - I had a brief, ill-guided idea to create a polygenic risk score from the rare variants of exome sequences that do not have externally determined effect sizes.  Rather than effect sizes, I would use computationally estimated severity scores (i.e. pathogenicity) along with annotation or some other tie-in to form the polygenic risk score.  The severity scores were not accurate enough to make anything useful.

**hypermobile** - Scripts that describe genome-wide association studies for various orthopedic conditions, with all data (and phenotypes) originating from the UK Biobank.

**jcCardio** - Investigation into rare, pathogenic variants and their association to cardiac rhythm conditions, namely atrial fibrillation.  We were able to define a set of variants that did signficantly associate to atrial fibrillation status, although not eventual adverse outcomes.

**why_wrong** - An attempt to determine why individuals with certain electronic health record defined disease labels (ICD-10, OPCS4) were linked to other information that would greatly suggest that they in fact do not have the disease.  No strong, workable solution was found to replace the disease labels, although many strong individual cases and general patterns were discovered.

**xwas** - Creation of polygenic risk scores with and without variants in the X-chromosome.  In this process the effect sizes that constitute the polygneic risk score were computed from a GWAS and the X-chromosome from an XWAS (with multiple options).  The sample size, after splitting the cohort, was too small to determine whether the X-added score was any better than the non X-added score.
