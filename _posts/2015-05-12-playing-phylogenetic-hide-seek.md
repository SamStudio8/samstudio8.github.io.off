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
early classifications appeared to focus on removing contradicting taxa
that were clouding "pure" definitions of what was truly plant and animal and
choosing to conviniently classify difficult organisms on what they are not, as
opposed to what they are[^5]. The term became a dumping ground for
unicellular protozoa that were not quite animals, protophytic algae that were
not quite plans and fungal organisms such as slime molds (it wouldn't be
until 1969 that Fungi would get a kingdom of their own).

After the invention of the first electron microscope in the early 20th century, confusion
was further compounded by both the appearance of a distinct cellular nucleus in some
unicellular organisms.

### The Protozoan Identity Crisis
Although Haeckel later added the now familiar term protozoa as a "sub-kingdom"
to his Protista kingdom to represent unicellular animals, another kingdom schema
emerged in 1938, proposed by Herbert Copeland in *The Kingdoms of Organisms*[^3].
Copeland moved the bacteria and algae out of Haeckel's Protista kingdom to create
a new kingdom: Monera.

Copeland would later re-name his kingdom of Protista to **Protoctista**
("first established beings") -- a term originally coined by John Hogg.
As Haeckel's models still persisted and the Protista kingdom still contained bacteria
it seemed "unfit" to continue using the same name in both models.
Selecting a new name, Copeland purposefully avoided the term Protozoa due
to both its confusing prior use as a kingdom, class and phylum, but also in
agreeance with a point made in Hogg's 1860 manuscript
*On the Distinctions of a Plant and an Animal and on a Fourth Kingdom of Nature*: 
 that the term "can alone include those that are admitted by all to be animals
 or 'zoa'"[^3]. *i.e.* The term is inappropriate if it is to be applied to non-animals.

### The Rise of Superkingdoms
In the 1960s, microbiologist Roger Stanier propagated the "fundamental division of life" observed
between the "prokaryotes" and "eukaryotes", originally noted decades prior by Edouard Chatton in 1925[^3].

In 1969, Robert Whittaker published his own five kingdom model (a revision of his earlier work
where he had in fact returned the bacteria to the protista kingdom, it would take a few more years of
research before concluding "this evolutionary divergence in cellular
structure [in bacteria] had to be accounted [for]"), but placing the kingdom Monera
beneath a new **Superkingdom of Bacteria** (placing the other kingdoms under
the new **Eukaryotic Superkingdom**) and introducing the kingdom of Fungi; a
wholly distinct kingdom from Plantae. Whittaker mainted that unicellularity was the most
important characteristic for membership of the kingdom of protista[^3].

Around the same time biologist Lynn Margulis began working on her own schema for organising life.
After several iterations, Margulis relied more on morphologic and structural observations over
Whittaker's unicellular criteria and allowed her Protista (later re-named to Protoctista in a nod
to Copeland) to contain eukaryotic
organisms that were "either unicellular **or multicellular** that are not plants, animals
or fungi"[Emphasis mine][^3].



### Sequencing
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
* Protozoa don't really exist and I need to be more specific about what I'm looking for.
* Even today, we still can't agree on what to call things and where they belong in a taxonomy, or even how best to present that taxonomy.


[^1]: Mark A. Ragan, *Trees and networks before and after Darwin*, 2009.

[^2]: Arshan Nasir, Kyung Mo Kim, and Gustavo Caetano-Anolles, *Giant viruses coexisted with the cellular ancestors and represent a distinct supergroup along with superkingdoms Archaea, Bacteria and Eukarya*, 2012.

[^3]: JM Scamardella, *Not plants or animals: a brief history of the origin of Kingdoms Protozoa, Protista and Protoctista*, 1999.

[^4]: AG Simpson and AJ Roger, *The real 'kingdoms' of eukaryotes*, 2004.

[^5]: <blockquote>"No mouth. No respiration. No entry."<br/><footer>— Kingdom <i>Animalia</i> Clubhouse Rules, 1860.</footer></blockquote>[^6]

[^6]: Based on paleontologist Richard Owen's 1858 description of plants and animals.
