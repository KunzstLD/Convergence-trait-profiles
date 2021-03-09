
---
<!-- bibliography: My Library.bib -->
<!-- output: word_document -->
---

- [Goal of the analysis:](#goal-of-the-analysis)
- [Methods](#methods)
  - [Statistical methods](#statistical-methods)
- [Results](#results)
  - [Trait profile groups across continents](#trait-profile-groups-across-continents)
  - [Importance of grouping features and traits for the clustering of macroinvertebrates](#importance-of-grouping-features-and-traits-for-the-clustering-of-macroinvertebrates)
- [Appendix](#appendix)


# Goal of the analysis:

1) Cluster taxa in the invertebrate trait databases according to their trait profiles to obtain groups of taxa with similar trait profiles (= trait profile groups)

2) Extract those traits that drive the clustering/are responsible for the grouping of the taxa in groups of similar trait profiles


# Methods 

Two analyses:

Analysis A) Comparison between Europe, North America, Australia, and New Zealand.

-  Aggregation to family level required because of varying taxonomic resolutions between trait databases

-   Subset of grouping features/traits that are present in all databases: feeding mode, locomotion, respiration, body size, body form, reproduction/oviposition, voltinism

Analysis B) Comparison between major climatic regions of Europe and North America:

-   Occurrence data for Europe based on Ecoregions

-   Occurrence data for North America from recently published CONUS trait database. Occurences are described for genera. Hence, aggregation of trait information to genus level required.


## Statistical methods

-   Data analysis:

    -   For Step 1) Hierarchical cluster analysis (taxa $\times$ trait matrix) to delineate trait profile groups (TPGs):

        -   Traits are expressed as affinities (based on fuzzy coding or binary coding)

        -   Normalised to a range between $[0-1]$

        -   Distance matrix based on a distance metric suitable for proportional data: Manly overlap index
        
        - GAP statistic to obtain the optimal number of groups   
        
        - Traits have been deemed as characteristic for their TPG, when more than 50 % of the families within a TPG express a trait with an affinity more than 0.5 (= defining traits). 

    -   For Step 2) Importance of grouping features and traits for the clustering of macroinvertebrates

        -    Contribution of the single variables to the global distance

        -    Random forest with groups as response; Permutation importance
        
             - Account for correlation among variables: RFE or Boruta? (Degenhardt paper) 
             - Search for trait interactions: iRF    

# Results


## Trait profile groups across continents

<!-- Add here: 
- Description of trait profile groups that occur in all regions
- In almost all regions
- Further general observations
- Include all heatmaps (can also go to the appendix) -->


The number of identified TPGs by hierarchical cluster analysis varied between continents. For Australia, 7 trait profile groups were delineated by the GAP statistic, for Europe 10, for North America 6, and for New Zealand 9. No TPG was identical on all continents, but three similar TPGs that shared certain traits could be delineated for all continents. In the following these TPGs that were similar in their defining traits across continents are denoted *converged group A* to *C* *(Kann man sicherlich später noch anders bennen)*. *Converged group A* consisted of gill-respirating crawlers, with a cylindrical body form that lay their eggs in an aquatic environment. *Converged group B* consisted of univoltine predators that lay their eggs in an aquatic environment. *Converged group C* consisted of small and cylindrical taxa, that lay their eggs in an aquatic environment

(*es gibt offensichtliche Ähnlichkeiten zwischen den TPGs, was bedeuten könnte das es weniger als drei ähnliche TPGs gibt. Ich denke insgesamt sind durch die GAP statstic zuviele Gruppen pro Kontinent abgeleitet wurden. Z.B. sind TPG3 und TPG4 in Australien sehr ähnlich und unterscheiden sich eig. nur dadurch das in einer Gruppe vermehrt Predatoren vorkommen (siehe die Heatmaps im Anhang). Vllt. wäre es besser noch einen anderen Weg zu wählen um Gruppen abzuleiten, oder dies a-posteriori festzulegen?*). 

![Heatmap displaying the traits of the *converged group A*. The defining traits for the *converged group A* are indicated by the orange rectangles. The left y-axis shows the families for each TPG. Continents and the number of the TPG delineated specific to the continent are displayed on the right y-axis. Other traits displayed are characteristic only for some of the displayed TPGs (e.g. TPG 6 in New Zealand was also characterised by univoltine, medium-sized predators).](Graphs/Heatmap_TPG1_across_all_continents.png)

![Heatmap displaying the traits of the *converged group A*. The defining traits for the *converged group A* are indicated by the orange rectangles. The left y-axis shows the families for each TPG. Continents and the number of the TPG delineated specific to the continent are displayed on the right y-axis. Other traits displayed are characteristic only for some of the displayed TPGs.](Graphs/Heatmap_TPG2_across_all_continents.png)

![Heatmap displaying the traits of the *converged group A*. The defining traits for the *converged group A* are indicated by the orange rectangles. The left y-axis shows the families for each TPG. Continents and the number of the TPG delineated specific to the continent are displayed on the right y-axis. Other traits displayed are characteristic only for some of the displayed TPGs](Graphs/Heatmap_TPG3_across_all_continents.png)


- Orders for the converged groups?


## Importance of grouping features and traits for the clustering of macroinvertebrates 

- Based on distance matrix: correlations between the (squared) distances obtained for each grouping feature and the global (squared) distances obtained by mixing all the grouping features (contributions of each grouping feature to the global distance)
  - For AUS: locomotion, oviposition, body form
  - EU: size, respiration, feeding mode
  - NOA: respiration, locomotion, body form
  - NZ: voltinism, size, respiration
- Respiration in three of four regions among the three most important grouping features (EU, NOA, NZ)


Table: Grouping features that contributed most to the global distance in the distance matrices in terms of correlation between the distances obtained for each grouping feature and the global distance obtained by mixing all the grouping features.

| Continent | Grouping feature | Correlation with <br />global distance |
| :-------- | :--------------- | -------------------------------------: |
| AUS       | locomotion       |                              0.5101570 |
| AUS       | oviposition      |                              0.4639916 |
| AUS       | body form        |                              0.4251654 |
| AUS       | size             |                              0.4239555 |
| AUS       | voltinism        |                              0.3905830 |
| AUS       | respiration      |                              0.3851332 |
| AUS       | feeding mode     |                              0.3500149 |
| EU        | size             |                              0.5166203 |
| EU        | respiration      |                              0.4881537 |
| EU        | feeding mode     |                              0.4523789 |
| EU        | voltinism        |                              0.4458472 |
| EU        | locomotion       |                              0.4083654 |
| EU        | oviposition      |                              0.2556404 |
| EU        | body form        |                              0.2482942 |
| NOA       | respiration      |                              0.4476170 |
| NOA       | locomotion       |                              0.4404900 |
| NOA       | body form        |                              0.4241282 |
| NOA       | size             |                              0.4043473 |
| NOA       | voltinism        |                              0.3938844 |
| NOA       | oviposition      |                              0.3656180 |
| NOA       | feeding mode     |                              0.3147170 |
| NZ        | voltinism        |                              0.5435340 |
| NZ        | size             |                              0.5417722 |
| NZ        | respiration      |                              0.4884413 |
| NZ        | feeding mode     |                              0.4422673 |
| NZ        | locomotion       |                              0.3815343 |
| NZ        | oviposition      |                              0.3164401 |
| NZ        | body form        |                              0.1688464 |


- Five most important traits initial run 
  - AUS: size_small, resp_gil, resp_pls_spi, feed_predator, ovip_aqu
  - EU: size_small, size_medium, feed_predator, locom_crawl, feed_herbivore
  - NOA: resp_gil, resp_teg, feed_herbivore, bf_cylindrical, size_large
  - NZ: feed_herbivore, size_small, locom_burrow, size_medium, volt_bi_multi
- Collinearity among variables, needs to be accounted for!
- *Train error looks like overfitting, test error relatively low -> indication for bad grouping by the clustering?*


Table: Five most important traits to distinguish TPGs for each continent. After the 
      initial run the highest ranking variable has been removed and a random forest has been fitted 
      again. This procedure has been repeated 4 times.

| Continent | Five most important traits                                                     | Run                    |
| :-------- | :----------------------------------------------------------------------------- | :--------------------- |
| AUS       | size_small, resp_gil, resp_pls_spi, feed_predator, ovip_aqu                    | initial                |
| AUS       | bf_cylindrical, size_medium   , bf_flattened  , ovip_ter      , resp_pls_spi   | rm_highest_ranking_var |
| AUS       | resp_pls_spi, ovip_aqu    , bf_flattened, size_medium , resp_gil               | rm_highest_ranking_var |
| AUS       | bf_cylindrical, size_medium   , resp_gil      , bf_flattened  , ovip_aqu       | rm_highest_ranking_var |
| AUS       | resp_pls_spi, ovip_aqu    , bf_flattened, size_medium , resp_gil               | rm_highest_ranking_var |
| EU        | size_small    , size_medium   , feed_predator , locom_crawl   , feed_herbivore | initial                |
| EU        | size_medium  , feed_predator, locom_swim   , locom_crawl  , resp_gil           | rm_highest_ranking_var |
| EU        | size_small    , locom_crawl   , feed_predator , feed_shredder , feed_herbivore | rm_highest_ranking_var |
| EU        | size_medium  , feed_predator, locom_swim   , locom_crawl  , resp_gil           | rm_highest_ranking_var |
| EU        | size_small    , locom_crawl   , feed_predator , feed_shredder , feed_herbivore | rm_highest_ranking_var |
| NOA       | resp_gil      , resp_teg      , feed_herbivore, bf_cylindrical, size_large     | initial                |
| NOA       | resp_teg      , feed_herbivore, locom_crawl   , bf_cylindrical, size_medium    | rm_highest_ranking_var |
| NOA       | resp_gil      , feed_herbivore, bf_cylindrical, size_large    , feed_predator  | rm_highest_ranking_var |
| NOA       | resp_teg      , feed_herbivore, locom_crawl   , bf_cylindrical, size_medium    | rm_highest_ranking_var |
| NOA       | resp_gil      , feed_herbivore, bf_cylindrical, size_large    , feed_predator  | rm_highest_ranking_var |
| NZ        | feed_herbivore, size_small    , locom_burrow  , size_medium   , volt_bi_multi  | initial                |
| NZ        | feed_predator, size_small   , locom_burrow , volt_bi_multi, size_medium        | rm_highest_ranking_var |
| NZ        | feed_herbivore, size_medium   , volt_bi_multi , size_small    , locom_swim     | rm_highest_ranking_var |
| NZ        | feed_predator, size_small   , locom_burrow , volt_bi_multi, size_medium        | rm_highest_ranking_var |
| NZ        | feed_herbivore, size_medium   , volt_bi_multi , size_small    , locom_swim     | rm_highest_ranking_var |

\pagebreak

# Appendix

![Heatmap displaying delineated trait profile groups (TPGs) for Australia. The grouping and order of the families on the y-axis was obtained from the dendrogram based on the hierarchical cluster analysis. The choice of the optimal number of groups was guided by the GAP statistic.](Graphs/Heatmap_tpgs_AUS.png)

![Heatmap displaying delineated trait profile groups (TPGs) for Europe. The grouping and order of the families on the y-axis was obtained from the dendrogram based on the hierarchical cluster analysis. The choice of the optimal number of groups was guided by the GAP statistic.](Graphs/Heatmap_tpgs_EU.png)

![Heatmap displaying delineated trait profile groups (TPGs) for North America. The grouping and order of the families on the y-axis was obtained from the dendrogram based on the hierarchical cluster analysis. The choice of the optimal number of groups was guided by the GAP statistic.](Graphs/Heatmap_tpgs_NOA.png)

![Heatmap displaying delineated trait profile groups (TPGs) for Europe. The grouping and order of the families on the y-axis was obtained from the dendrogram based on the hierarchical cluster analysis. The choice of the optimal number of groups was guided by the GAP statistic.](Graphs/Heatmap_tpgs_NZ.png)