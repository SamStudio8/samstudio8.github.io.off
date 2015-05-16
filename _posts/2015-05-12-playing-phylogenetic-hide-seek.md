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

<blockquote>[...] single celled micro-organisms that feed from their direct surroundings
and have the capacity for controlled movement with a tendency to thrive in moist
environments [...]</blockquote>

Yet it seems the answer to "What are protozoa?" boils down to who you ask, or rather,
whose interpretation of the taxonomic system you ask.


## A Very Brief History of Modern Taxonomy
In the 18th century, Carl Linnaeus published works on the biological classification of organisms based
on their structural appearances, establishing the three kingdoms of *Animale* (Animal), *Vegetabile* (Vegetable)
and *Lapideum* (Mineral) which lend their names to a guessing game of the same name.
Linnaeus' work, particularly on naming strategies for organisms served as a foundation to modern taxonomy.

However, the class system set out by Linnaeus was designed for the identification of already catalogued species.
Following Charles Darwin's revolutionary *On the Origin of Species* a century later, tree representations of
species sharing common descent because increasingly common (although not the first). Throughout the late
19th and early 20th century scientists attempted to group available fossil records by structual affinity
to link ancestry of common species through the ages. It was here we drew the comparison between
contemporary bird species and the dinosaurs.


Thomas Cavalier-Smith proposed in 1998 an oft-used "six kingdom"

With the development of chain-termination DNA sequencing in 1977 by Frederick Sanger, the field
of molecular genetics was born. Taxonomy could now be based on differences expressed in the genetic
sequences of organisms (particularly in highly conserved subsequences, such as those responsible
for constructing ribosomal RNA molecules), rather than by subjective observations of appearance and function.

However, towards the end of the 80s, as the fields of molecular biology and
In 1989, Michael Sleigh, in the second edition of *Protozoa and other Protists* remarked:

<blockquote>
The position is now clearer, and there is much support for the view that eukaryotes are best divided
into four kingdoms: Animalia or multicellular animals [...], Plantae or green land plants [...], Fungi
and Protista, comprising eukaryoute groups formerly classified as algae, protozoa and flagellate fungi.
<footer>— Michael Sleigh, *Protozoa and other Protists* (2nd ed.), 1989.</footer>
</blockquote>

Indeed, in the first edition of this book, *The Biology of Protozoa*, published in 1973 -
whose preface describes the "flourishing"
state of research in to Protozoa; "these organisms provide excellent subjects for studies on general
biological phenomena at the cellular level", primarily due to the ease of which large populations
of similar cells at the same age can be obtained, uncontaminated. 

*Evolution: Revisiting the Root of the Eukaryote Tree*

...currently the domain model holds no virii (but they have to go somewhere...) as they
are classed at non-cellular... though ...there is support for adding a fourth domain of life for virii[^2]

* * *

# tl;dr
* Even today, we still can't agree on what to call things and where they belong in a taxonomy, or even how best to present that taxonomy.
* 

[^1]: Mark A. Ragan, *Trees and networks before and after Darwin*, 2009.
[^2]: Arshan Nasir, Kyung Mo Kim, and Gustavo Caetano-Anolles, *Giant viruses coexisted with the cellular ancestors and represent a distinct supergroup along with superkingdoms Archaea, Bacteria and Eukarya*, 2012.
