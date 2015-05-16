---
layout: post
title: "Playing Phylogenetic Hide and Seek with Protozoa [WIP]"
---

Amanda suggested that alongside archaeal, bacterial and fungal associated hydrolases,
we should also look at [**protozoans**](http://en.wikipedia.org/wiki/Protozoa). 
No problem, I'll just get the taxonomy ID for protozoa and extract another database
from UniProtKB [as before]({% post_url 2015-04-24-trembling %}). Simple! Or so I thought...

Classification of protozoa appears to be less clear than I had realised. UniProtKB
[lists only three taxonomy entries](http://www.uniprot.org/taxonomy/?query=protozoa&sort=score)
for the term:

* uncultured rumen protozoa
* uncultured Canadian Arcott wether rumen protozoa
* uncultured protist

But none of these are really what I'm looking for. UniProtKB uses these psuedo-species
as a catch-all for environmental samples that don't really fit elsewhere in the database.
but I want a high level [taxonomic rank](http://en.wikipedia.org/wiki/Taxonomic_rank) like a kingdom
that encompasses all of the organisms of interest. So what are the organisms of interest?
In an early introduction to my project, I described the protozoa as:

<blockquote>[...] single celled micro-organisms that feed from their direct surroundings and have the capacity for controlled movement with a tendency to thrive in moist environments [...]</blockquote>

Yet it seems the answer to "What are protozoa?" boils down to who you ask, or rather,
whose interpretation of the taxonomic system you ask.

## A Very Brief History of Modern Taxonomy
In the 18th century, Carl Linnaeus published works on the biological classification of organisms based
on their structural appearances, establishing the three kingdoms of *Animale* (Animal), *Vegetabile* (Vegetable)
and *Lapideum* (Mineral) which lend themselves to a guessing game of the same name.
Linnaeus' work, particularly on [naming strategies for organisms](http://en.wikipedia.org/wiki/Binomial_nomenclature) served as a foundation to modern taxonomy.

However, the class system set out by Linnaeus was designed for identification rather than heritage.
Following Charles Darwin's revolutionary *On the Origin of Species* a century later in 1859, tree
representations of species sharing common descent became increasingly common
(although his tree diagrams were by no means the first) and as evolutionary theory
developed and took hold, scientists attempted to group available fossil
records by structual affinity to link ancestry of common species through the ages.
It was here we drew the comparison between contemporary bird species and the dinosaurs.

Meanwhile, micro-organisms exhibiting a variety of characteristics with affinity to both known
plants and animals were confusing the topic of taxonomy further...

### All Hail the Mighty Kingdom of Protista
Before Darwin, in 1820, the first use of **Protozoa** appeared in literature from Georg Goldfuss as
a class for "first, or early animals" inside the established kingdom of *Animalia*. Other scientists
believed that the protozoa were in fact a phylum, or even a kingdom of their own[^3].

In 1866, Ernst Haeckel proposed removal of the "mineral" kingdom (as Haeckel
did not recognise it as a kingdom of life) and the addition of the kingdom of **Protista** to
contain "doubtful organisms of the lowest rank which display no decided affinities nearer to
one side [animals] than to the other [plants]" or organisms possessing "animal and vegetable
characters united and mixed"[^3]. Haeckel proposed this as a "kingdom of primitive forms"
which included the "Monera" (bacteria) as members of this kingdom too.

It's clear how such kingdoms are later described in 21st century literature
as "a grab-bag for all eukaryotes that are not animals, plants or fungi"[^4], these
early classifications appeared to focus on removing confusing taxa
that were clouding "pure" definitions of what was truly plant and animal and
choosing to conviniently classify difficult organisms on what they are not, as
opposed to what they are[^5].


Throughout this period in the late 19th and early 20th century, various interpretations on how
species should be organised were proposed. Most schemes used the term "kingdom"


With the development of chain-termination DNA sequencing in 1977 by Frederick Sanger, the field
of molecular genetics was born. Taxonomy could now be based on differences expressed in the genetic
sequences of organisms (particularly in highly conserved subsequences, such as those responsible
for constructing ribosomal RNA molecules), rather than by subjective observations of morphology and function.

Different sequence database began to create their own taxonomic classification structure


In 1989, Michael Sleigh, in the second edition of *Protozoa and other Protists* remarked:

<blockquote>
The position is now clearer, and there is much support for the view that eukaryotes are best divided
into four kingdoms: Animalia or multicellular animals [...], Plantae or green land plants [...], Fungi
and Protista, comprising eukaryoute groups formerly classified as algae, protozoa and flagellate fungi.
<br/>
<footer>— Michael Sleigh, <i>Protozoa and other Protists</i> (2nd ed.), 1989.</footer>
</blockquote>

Indeed, in the first edition of this book, *The Biology of Protozoa*, published in 1973 -
whose preface describes the "flourishing"
state of research in to Protozoa; "these organisms provide excellent subjects for studies on general
biological phenomena at the cellular level", primarily due to the ease of which large populations
of similar cells at the same age can be obtained, uncontaminated. 

## Domains and Kingdoms

Thomas Cavalier-Smith proposed in 1998 an oft-used "six kingdom"
*Evolution: Revisiting the Root of the Eukaryote Tree*

...currently the domain model holds no virii (but they have to go somewhere...) as they
are classed at non-cellular... though ...there is support for adding a fourth domain of life for virii[^2]...


## UniProtKB
So how does this all relate to UniProt? The help documentation offloads the responsibility to the NCBI:

<blockquote>The taxonomy database that is maintained by the UniProt group is based on the <a href="http://www.ncbi.nlm.nih.gov/taxonomy">NCBI taxonomy database</a>, which is supplemented with data specific to the UniProt Knowledgebase (UniProtKB). While the NCBI taxonomy is updated daily to be in sync with GenBank/EMBL-Bank/DDBJ [...]</br><footer>— <a href="http://www.uniprot.org/help/taxonomy">UniProt Help: Taxonomy</a></footer></blockquote>

Chapter 4 of the [*NCBI Handbook*](http://www.ncbi.nlm.nih.gov/books/NBK21100/), *The Taxonomy Project*
explains the phylogenetic taxonomy


* * *

# tl;dr
* Even today, we still can't agree on what to call things and where they belong in a taxonomy, or even how best to present that taxonomy.

[^1]: Mark A. Ragan, *Trees and networks before and after Darwin*, 2009.

[^2]: Arshan Nasir, Kyung Mo Kim, and Gustavo Caetano-Anolles, *Giant viruses coexisted with the cellular ancestors and represent a distinct supergroup along with superkingdoms Archaea, Bacteria and Eukarya*, 2012.

[^3]: JM Scamardella, *Not plants or animals: a brief history of the origin of Kingdoms Protozoa, Protista and Protoctista*, 1999.

[^4]: AG Simpson and AJ Roger, *The real 'kingdoms' of eukaryotes*, 2004.

[^5]: <blockquote>"No mouth. No respiration. No entry."<br/><footer>— Kingdom <i>Animalia</i> Clubhouse Rules, 1860.</footer></blockquote>[^6]

[^6]: Based on paleontologist Richard Owen's 1858 description of plants and animals.
