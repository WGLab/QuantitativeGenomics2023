# Introduction
This repository contains the course materials for hands-on exercise of the "Variants Annotation and Phenotype Analysis" session in the Quantitative Genomics workshop at Columbia University on June 14-15, 2023.

# Preparation of computing environment
Both the lectures and hands-on exercises will be taught via Zoom video conference online. To ensure the cross-platform compatibility, we will only use Rstudio Cloud to implement tools that are developed in Perl and Python.

Students need to go to [Posit Cloud](https://posit.cloud/) to create your own account before the lectures: 
1. click 'Sign up' on the top-right corner
![image](https://github.com/WGLab/QuantitativeGenomics2023/assets/11565618/6cbaa160-c984-46eb-8960-8783d3a84a37)

2. select "Free" and click "Sign up" at the right panel. 
![image](https://github.com/WGLab/QuantitativeGenomics2023/assets/11565618/e681015d-7d7c-4bce-a766-936f78b56d63)

3. After you input your information, click "Sign up". 
<img src="https://user-images.githubusercontent.com/16017780/122824633-f0248600-d2ae-11eb-9bd1-a72b8f03383b.png" width="200" height="300">

If you logout after signing up, you can go to [RStudio Cloud](https://rstudio.cloud/) and login with your username and password. After you are added to "2023 Quantitative Genomics" project, you can click "Annotation_Assignment" to start the exercise. The direct link is [here](https://rstudio.cloud/spaces/232505/content/4126775). 

![image](https://user-images.githubusercontent.com/5926328/173241575-fa2b97b0-072e-4dc4-8876-d4ed33197ca5.png)

When the student opens an assignment project, RStudio Cloud will automatically make a copy for them. In other words, each student will work on your own copy of the virtual machine without influencing each other's files.

By default, you will be in "Console" as shown below where you can run R code.

![image](https://user-images.githubusercontent.com/5926328/173241769-fcf70153-9f11-4075-97a1-f44384c09d80.png)

You can also click "Terminal" so that you can run linux commands.

![image](https://user-images.githubusercontent.com/5926328/173241835-72cd8679-e36b-496f-ad91-98d2e35fc1af.png)

During the exercise, you can feel free to swith to Terminal when you need to run linux commands and to Console when you need to run R code.

# Basic Linux commands

After login to Rstudio Cloud, you already have terminal access by clicking 'terminal' button and can run all basic Linux commands. If you are not familiar with basic Linux commands, you can follow one simple tutorial to learn several commands: https://maker.pro/linux/tutorial/basic-linux-commands-for-beginners. Some of the  commands that we will use in the exercise include `ls`, `pwd`, `cd`, `mv` and `wget`.

# Functional annotation of genetic variants

### 1. Install ANNOVAR

Typically you will go to the [ANNOVAR website](http://annovar.openbioinformatics.org), fill in a registration form, and download the package there. For this exercise, we already prepared a folder in the cloud server that contains a "compact" version of ANNOVAR and necessary library files, to make it easier for users. The folder is called `genomics_exercise` and you can simply use `cd genomics_exercise` to go to this folder). To confirm which directory you are in, you can type `pwd`. You will see the example below.

```
/cloud/project$ cd genomics_exercise/
/cloud/project/genomics_exercise$ pwd
/cloud/project/genomics_exercise
```

Please note that this is a sub-sampled version of ANNOVAR for the purpose of th exercise today (to reduce file size significantly to <200Mb), and it should not be used in any other real data analysis. You can view all the files in `exercise1` directory:

![image](https://user-images.githubusercontent.com/5926328/173242142-864ba9a2-5ac2-4d9f-94f8-a3c4a3d61d52.png)


### 2. Run ANNOVAR on a small VCF file

Type `cd exercise1` to enter the `exercise1` directory. The sub-folder `humandb` folder already contains several annotation databases for human genome that we will use in our exercise. (Note that users can find more annotation databases [here](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/user-guide/download/#-for-filter-based-annotation).

```
perl table_annovar.pl example/ex2.vcf humandb/ -buildver hg19 -out myanno -remove -protocol refGeneWithVer,cytoBand,gnomad211_exome -operation g,r,f -nastring . -vcfinput -polish
```

After that, you will find the result files `myanno.hg19_multianno.txt` and `myanno.hg19_multianno.vcf`.

The `-protocol` and the `-operation` arguments are the most important one to understand here. Protocol lists a set of tasks that you want to use for annotation, composed of a list of comma-separated keywords. Operation specifies the type of annotations that each step in the protocol uses, composed of a list of comma-separated letters such as `g`, `r` and `f`, corresponding to gene-based, region-based and filter-based annotations, respectively.

In the command above, we specified a protocol with three annotation tasks, including RefGene annotation (a gene-based annotation), cytogenetic band annotation (a region-based annotation), and allele frequency in gnoMAD version 2.1.1 (a filter-based annotation). These three types of annotation are indicated as `g`, `r` and `f` in the -operation argument, respectively. We also specify that dots ('.') be used to indicate lack of information, so variants that are not observed in the database will have '.' as the annotation. The `-vcfinput` argument specifies that the input file is in VCF format. We also specify that the output file names should have "myanno" as the prefix, to make it easier to check output files.

If you have Excel installed, you can open the `hg19_multianno.txt` file by Excel and examine the various tab-delimited fields. Since we are using Rstudio now, we could also easily import the tab-delimited file into Rstudio to examine the various rows and columns. We will do this later to filter and prioritize disease-relevant variants.

### 3. Run ANNOVAR on a human exome

Next, we want to download a VCF file and then run ANNOVAR on this file.

In the `exercise1` folder:

```
wget http://molecularcasestudies.cshlp.org/content/suppl/2016/10/11/mcs.a001131.DC1/Supp_File_2_KBG_family_Utah_VCF_files.zip -O Supp_File_2_KBG_family_Utah_VCF_files.zip
```

To give some background information, this is a zip file as supplementary material of a published paper on exome sequencing of a family with undiagnosed genetic diseases. Through analysis of the exome data, the proband was confirmed to have KBG syndrome, a disease caused by loss of function mutations in ANKRD11 gene. There are several VCF files contained in the zip file, including those for parents, silings and the proband. We will only analyze proband in this exercise, but if you are interested, you may want to check whether this is a de novo mutation by analyzing parental genomes.

Next unzip the ZIP file (`unzip`):

```
unzip Supp_File_2_KBG_family_Utah_VCF_files.zip
mv './File 2_KBG family Utah_VCF files/' VCF_files
```

Note that in the commands above, we rename the directory to "VCF_files" just to make it easier to read and operate.

Run ANNOVAR on the VCF file:
```
perl table_annovar.pl VCF_files/proband.vcf humandb/ -buildver hg19 -out proband.annovar -remove -protocol refGeneWithVer,gnomad211_exome -operation g,f -nastring . -vcfinput
```

The `proband.annovar.hg19_multianno.txt` file contains annotations for this exome.

Now compare this command from the previous table_annovar.pl command, we can see that this type we only request two annotation tasks (refGeneWithVer and gnomad211_exome), which are gene-based and filter-based annotations, respectively.


### 4. Results visualization

We used "terminal" for the commands above to generate results. To do further analysis in R, we can now switch to the "console" in Rstudio.

In the `exercise1` folder, run `pwd` to obtain the working directory and copy this path. 

Paste this path to R console or R studio interface to setup the working directory using function `setwd()`. For example:

```
setwd("/cloud/project/genomics_exercise/exercise1")
```

Check variant distribution across chromesomes:
```{r}
#load data
res <- read.table("proband.annovar.hg19_multianno.txt", fill=T, header=T, sep="\t", na.strings = ".")

#visualize variant frequency
attach(mtcars)
par(mar=c(5.1, 4.1, 4.1, 2.1),mfrow=c(1,1))
table <- table(res$Chr)
chrlist <- paste0("chr",1:22)
chrlist <- c(chrlist, "chrX")
barplot(table[chrlist], ylab="Number of Variant", las=2)
````

At this point, you should see a barplot similar to the one below:

![image](https://user-images.githubusercontent.com/11565618/123002115-87590e80-d37f-11eb-8b90-bf2901e3a8c8.png)

Check allele frequency distribution between non-synonymous, synonymous and intronic variants:
```{r}
#visualize allele frequency
par(mar=c(2, 4, 1, 1),mfrow=c(3,1))
af_syn<-res[res$ExonicFunc.refGeneWithVer=="synonymous SNV",]$AF
hist(as.numeric(af_syn),  ylab="Frequency", xlab="Allele frequency", main="Synonymous SNV", breaks=seq(0,1,.001))

af_nonsyn<-res[res$ExonicFunc.refGeneWithVer=="nonsynonymous SNV",]$AF
hist(as.numeric(af_nonsyn), ylab="Frequency", xlab="Allele frequency", main="nonSynonymous SNV",breaks=seq(0,1,.001))

af_intron<-res[res$Func.refGeneWithVer =="intronic",]$AF
hist(as.numeric(af_intron), ylab="Frequency", xlab="Allele frequency", main="Intronic", breaks=seq(0,1,.001))
```

You should see a figure similar to the one below:

![image](https://user-images.githubusercontent.com/11565618/123002219-a9eb2780-d37f-11eb-9452-01cf20f5eedc.png)



Check allele frequency distribution across race:
```{r}
#visualize allele frequency across race
par(mar = c(12, 5, 2, 2),mfrow=c(1,1)) 
res_sub <- res[res$Gene.refGeneWithVer=="NRXN2",]
variantname <- res_sub$Otherinfo6[1]
res_sub <- as.matrix(res_sub[1,16:23])
colnames(res_sub) <- c("African/African-American", "South Asian", 
                       "Latino/Admixed American", "East Asian", "Non-Finnish European", "Finnish",
                       "Ashkenazi Jewish", "Other")
barplot(res_sub, ylab="Allele frequency", main=variantname, las=2)
```

You should see a figure similar to the one below:

![image](https://user-images.githubusercontent.com/11565618/123002260-ba030700-d37f-11eb-811b-54eb98fb0c73.png)

Feel free to use the table browser to examine the various rows and columns in this table and perform additional summary statistics. Later on, after the exercise 2 (phenotype analysis) below, we will go back to this table to show how combined analysis of genotype data and phenotype data can facilitate genetic diagnosis of rare diseases.


# Phenotype-driven prioritization of human disease genes (Phen2Gene, ClinPhen, AMELIE, etc)

Phen2Gene is a phenotype-based gene prioritization tool from HPO IDs or clinical notes on patients. In the next exercise, we will first use Phen2Gene to prioritize genes based on clinical phenotypes of patients with Mendelian diseases.

We will do the exercise in a directory called `exercise2`. In the terminal, assuming that you are currently in the `exercise1` directory (you can use `pwd` to check this), you can use command `cd ../exercise2` to go into the exercise2 directory.

There are three ways to run Phen2Gene: download and run it locally (need a few GB of space), using 

and using Phen2Gene website.

The benefit of running Phen2Gene is if you do not have any idea of a candidate gene for your disease, you can use it in one of three scenarios:

1. Ideally, you have a list of physician-curated HPO terms describing a patient phenotype and a list of potential candidate disease variants that overlap gene space and you want to narrow down the list of variants by prioritizing candidate disease genes, often in tandem with variant prioritization software, which cannot as of yet score STR expansions or SVs unlike Phen2Gene which is variant agnostic.
2. You do not have variants, but you have HPO terms and would like to get some candidate genes for your disease that you may want to target sequence, as it is much cheaper than whole-genome or whole-exome sequencing.
3. If you have clinical notes, you can use tools like EHR-Phenolyzer or Doc2HPO for processing clinical notes into HPO terms using natural language processing (NLP) techniques, then apply scenario 1 or 2 as relevant.

Other tools listed below (ClinPhen, AMELIE, etc) require a gene list, and Phen2Gene does not require any variant data or prior gene lists to provide high quality candidate genes.  One would most likely use Phen2Gene for obtaining the genetic etiology of an undiagnosed rare disease with no obvious genetic cause.


### 1. Using Phen2Gene API
1. Go to Terminal, make sure you are in the `exercise2` directory first, and run `curl -i -H "Accept: application/json" -H "Content-Type: application/json" "https://phen2gene.wglab.org/api?HPO_list=HP:0000455;HP:0000574;HP:0030084;HP:0012471;HP:0000239;HP:0001572;HP:0000960;HP:0001250;HP:0000322;HP:0001831" | tail -n 1 > output.txt`
where you generate JSON output in `output.txt`
However, since the `output.txt` file is in JSON format, it is not very intuitive to view the content of the file. Instead, we will use the table browser in Rstudio to view a formatted version of the JSON file.
2. Go To Console, remember that we are probably in the `exercise1` directory, so we should first set `exercise2` as the working directory.
```
setwd("../exercise2")
```
3. Then we can run
```
# install JSON in R
install.packages("rjson")
# Load JSON in R
library("rjson")
# Read JSON results
result <- fromJSON(file = "output.txt")
# Convert them to array and name columns.
marray <- t(array( unlist(result$results), dim=c(5, length(result$results)) ) );
colnames(marray) <- names(result$results[[1]]);
# View the results in 2-D array. The second column is the rank of genes.
View (marray);
# Get top 100 genes based on Phen2Gene score
top100<-matrix(array(unlist(result$results)), 5)[1,][1:100]
write.table(top100, file="phen2gene_top100_genes", quote=F, col.names = F, row.names = F)
```
![image](https://user-images.githubusercontent.com/5926328/173119713-fd122a2a-8b4d-42d8-9618-af9b06ad5e8e.png)

You can see that the top ranked genes are VPS13B, ARID1B, etc.

### 2. Running Phen2Gene locally in Terminal

Under the `exercise2` directory, we can run `python /cloud/project/usr/Phen2Gene/phen2gene.py -f hpo_list.txt`
where we generate Phen2Gene output file in `out/output_file.associated_gene_list`.

The `hpo_list.txt` file contains a list of HPO terms for a patient (one per line). If you want to see what is in the file, you can use `cat hpo_list.txt` to see the list of terms. The output file is written to `out/output_file.associated_gene_list` by default.

Use command `more out/output_file.associated_gene_list` to check the predicted genes.

![image](https://user-images.githubusercontent.com/11565618/173251477-9995c724-a637-4a81-8582-31ac593ec401.png)

You can see that the top ranked genes are VPS13B, ARID1B, etc. Type `q` to quit the mode.

### 3. Using the Phen2Gene Website to assess the case with KBG syndrome

Go to https://phen2gene.wglab.org.  Click on the tab `Patient notes`:

![image1](https://user-images.githubusercontent.com/6568964/84083556-df5ce700-a9af-11ea-87c5-d02cc742b8c4.png)

and paste this paragraph of clinical notes in the text box:

```
The proband had an abnormally large fontanelle, which resolved without treatment. The proband does not appear to have a sacral dimple. Other than the presence of flat arches, there are no obvious signs of foot abnormalities. The proband does not look like his other siblings, although there was a resemblance between him and his sister when they were the same age. Features of proband’s face can be seen, including bushy eyebrows, broad nasal tip, short philtrum, full lips and cupid bow of upper lip. Video of proband’s behavior. Proband is non-verbal, and hyperactive. He repetitively spins his toy. While playing, he gets up from his chair, walks a few steps, stomps his feet, and sits back down. Additional interview in August 2015, after the parents received the proband’s diagnosis of KBG syndrome. The proband develops new mannerisms every four to six months, the most recent being short, hard breaths through the nose and head turning. The proband has had a substantial decrease in the number of seizures after starting an Epidiolex (cannabidiol) treatment (70-80% decrease as described by the parents). The frequency of seizures increased after the proband fell and fractured his jaw.  The mother describes the proband’s macrodontia. Although the mother and several siblings have large teeth, the macrodontia in the proband does not appear in any other member of the family.  The proband’s features are compared to other characteristics usually found in other KBG patients. Unlike most KBG patients, the proband has full lips. Like most KBG patients, the proband has curved pinkies (diagnosed as clinodactyly), which are often found in KBG patients.  Although the proband has relatively short toes, this trait may have been inherited from the father. The proband also has curved toenails, which commonly appear in autistic children.
```

![image2](https://user-images.githubusercontent.com/6568964/84083597-f26fb700-a9af-11ea-8fff-ead80c87d479.png)

Then click Submit.

![image3](https://user-images.githubusercontent.com/6568964/84083616-fac7f200-a9af-11ea-97cd-68585539d7fe.png)

You should see that ANKRD11 is in the top 3.  The reason it is not 2 as in the previous example is that we manually curated HPO terms from the patient notes for each patient in our Phen2Gene paper and this is raw notes with no curation done and all HPO terms are being automatically extracted by Doc2HPO.  Despite this, the discrepancy in our results is minimal.

![image4](https://user-images.githubusercontent.com/6568964/84083649-0fa48580-a9b0-11ea-81ed-7d999f40eade.png)

Alternatively, you can also submit the HPO terms from `example/ANKRD11.txt` manually using the tab `HPO IDs` (they are already the default terms in the website so you can just click `Submit` on that tab).

![image5](https://user-images.githubusercontent.com/6568964/84083556-df5ce700-a9af-11ea-87c5-d02cc742b8c4.png)

And you should see that ANKRD11 is still number 2.

![image6](https://user-images.githubusercontent.com/6568964/84211809-3afba300-aa8a-11ea-8674-89518b0f8576.png)

### 4. Run ClinPhen

ClinPhen is another tool to analyze clinical notes. To run this tool, 
1. Go to Terminal and run the commands below to install it.
```

export PATH="$HOME/.local/bin":$PATH

pip install clinphen
```

2. Download package nltk in Python. In terminal, type `python`, then type
```
import nltk
nltk.download('omw-1.4')
```
After completing installation, type `exit()` to quit Python.

3. Run ClinPhen
```
clinphen m_notes.txt
```

The expected results are below:

```
/cloud/project/genomics_exercise/exercise2$ clinphen m_notes.txt
HPO ID  Phenotype name  No. occurrences Earliness (lower = earlier)     Example sentence
HP:0012471      Thick vermilion border  2       12      unlike most kbg patients the proband has full lips 
HP:0000574      Thick eyebrow   1       9       features of proband s face can be seen including bushy eyebrows broad nasal tip short philtrum full lips and cupid bow of upper lip 
HP:0000455      Broad nasal tip 1       10      features of proband s face can be seen including bushy eyebrows broad nasal tip short philtrum full lips and cupid bow of upper lip 
HP:0000322      Short philtrum  1       11      features of proband s face can be seen including bushy eyebrows broad nasal tip short philtrum full lips and cupid bow of upper lip 
HP:0002263      Exaggerated cupid's bow 1       13      features of proband s face can be seen including bushy eyebrows broad nasal tip short philtrum full lips and cupid bow of upper lip 
HP:0001250      Seizures        1       32      the frequency of seizures increased after the proband fell and fractured his jaw 
HP:0030084      Clinodactyly    1       42      like most kbg patients the proband has curved pinkies diagnosed as clinodactyly which are often found in kbg patients 
```

You will get the results from ClinPhen. You can also type `clinphen --help` if you are interested in more options.

To find candidate genes for this set of HPO terms, you can input `HP:0012471;HP:0000574;HP:0000455;HP:0000322;HP:0002263;HP:0001250;HP:0030084` into Phen2Gene web server.

### 5. Running AMELIE

[AMELIE](https://amelie.stanford.edu/) is yet another tool for analyzing clinical notes, however it requires a candidate gene list.  We put the real causal gene in a algorithm-allowed maximum 1000 gene list of random exome capture genes.

First copy [this list of HPO terms](https://raw.githubusercontent.com/WGLab/Phen2Gene/master/example/HPO_sample.txt) to your clipboard. Then go to the [AMELIE website](https://amelie.stanford.edu/) and click "Submit Case."

![image](https://user-images.githubusercontent.com/6568964/122857269-3812cf80-d2e6-11eb-9649-4a30f0d8d94d.png)

Now, paste the HPO term list into "Case phenotypes."

![image](https://user-images.githubusercontent.com/6568964/122857504-93dd5880-d2e6-11eb-94b8-d4ce749df014.png)

Then, copy [this list of genes](https://raw.githubusercontent.com/WGLab/Phen2Gene/master/example/1000genetest.txt) to your clipboard.  

Next, paste the gene list into "Case genotype" under the tab "Gene List" then click Submit.

![image](https://user-images.githubusercontent.com/6568964/122857542-a2c40b00-d2e6-11eb-82ca-ac10ea712772.png)

You'll get a screen that looks like this:

![image](https://user-images.githubusercontent.com/6568964/122857711-ea4a9700-d2e6-11eb-8c03-761257d6ff8b.png)

And a wheel that tells you to wait a _while_.

After the wheel is done spinning your results should look like the following:

![image](https://user-images.githubusercontent.com/6568964/122857763-051d0b80-d2e7-11eb-8e14-2941636be222.png)

The number 1 gene, TAF1 is the correct causal gene.  But Phen2Gene gets the same result.

### 6. Running PhenCards

Go to https://phencards.org/.  Click on the tab `Patient notes`:

![image](https://user-images.githubusercontent.com/11565618/172928326-9406c762-9887-46d3-80cf-b8098c88bb27.JPG)

and paste this paragraph of clinical notes in the text box:

```
The proband had an abnormally large fontanelle, which resolved without treatment. The proband does not appear to have a sacral dimple. Other than the presence of flat arches, there are no obvious signs of foot abnormalities. The proband does not look like his other siblings, although there was a resemblance between him and his sister when they were the same age. Features of proband’s face can be seen, including bushy eyebrows, broad nasal tip, short philtrum, full lips and cupid bow of upper lip. Video of proband’s behavior. Proband is non-verbal, and hyperactive. He repetitively spins his toy. While playing, he gets up from his chair, walks a few steps, stomps his feet, and sits back down. Additional interview in August 2015, after the parents received the proband’s diagnosis of KBG syndrome. The proband develops new mannerisms every four to six months, the most recent being short, hard breaths through the nose and head turning. The proband has had a substantial decrease in the number of seizures after starting an Epidiolex (cannabidiol) treatment (70-80% decrease as described by the parents). The frequency of seizures increased after the proband fell and fractured his jaw.  The mother describes the proband’s macrodontia. Although the mother and several siblings have large teeth, the macrodontia in the proband does not appear in any other member of the family.  The proband’s features are compared to other characteristics usually found in other KBG patients. Unlike most KBG patients, the proband has full lips. Like most KBG patients, the proband has curved pinkies (diagnosed as clinodactyly), which are often found in KBG patients.  Although the proband has relatively short toes, this trait may have been inherited from the father. The proband also has curved toenails, which commonly appear in autistic children.
```

![image](https://user-images.githubusercontent.com/11565618/172928866-dfde622c-3a00-42dc-9eab-0dba303ff28f.JPG)

Then click Search.

![image](https://user-images.githubusercontent.com/11565618/172929187-2d6ea4d2-a333-4ee6-8914-5b6b5db0b9d7.JPG)

This is raw notes with no curation done and all HPO terms are being automatically extracted by Doc2HPO. You should see that KBG syndrome is among the top disease list.

![image](https://user-images.githubusercontent.com/5926328/173098943-9e24c858-318b-45f6-9ca9-ea3d1ea1c0aa.png)

Alternatively, you can also submit the following list of HPO terms manually using the PhenCards web server to find the information of HPO terms. (Do not enter these terms into PhenCards which only find cross-reference for each one of these terms.)

```
HP:0000455; HP:0000574; HP:0030084; HP:0012471; HP:0000239; HP:0001572; HP:0000960; HP:0001250; HP:0000322; HP:0001831
```

The expected results are shown below in Phen2Gene server:

![image](https://user-images.githubusercontent.com/5926328/173100269-d316d071-c96e-46d4-b28e-9c479ceb7e9c.png)

### 7. Filtering ANNOVAR annotation results based on Phen2Gene output
Now that we performed variant annotation in exercise1, and performed phenotype-based gene prioritization in exercise2, we want to see how to infer and prioritize candidate genes for the disease.

Go to Console, we can run
```
# Load top 100 genes generated by Phen2Gene
top100 <- read.table("out/output_file.associated_gene_list", header=T)
top100 <- top100$Gene[1:100]

# Load results of ANNOVAR
res <- read.table("../exercise1/proband.annovar.hg19_multianno.txt", fill=T, header=T, sep="\t", na.strings = ".")

# Filtering ANNOVAR annotations based on Phen2Gene output
res$Gene.refGeneWithVer
res <- res[res$Gene.refGeneWithVer %in% top100,]
res <- res[res$Func.refGeneWithVer=="exonic",]
res <- res[res$ExonicFunc.refGeneWithVer!="synonymous SNV",]
res <- res[res$AF_popmax<.005 | is.na(res$AF_popmax),]
View(res)
```

The analytical logic in the above command is that: first we take only exonic variants in the top 100 genes predicted by phenotype-driven gene prioritization software. Next, we remove synonymous SNVs from the list. Next, we keep only variants that are not observed in the gnomAD database, or observed but with allele frequency lower than 0.005.

2. You can see that gene ANKRD11 is in the following list

![image](https://user-images.githubusercontent.com/11565618/173124914-4cc9f913-c265-4153-ba8c-16d7ec66f716.png)

### 8. Running Phenomizer

[Phenomizer](https://compbio.charite.de/phenomizer/) is another a web-based application for clinical diagnostics in human genetics using semantic similarity searches in ontologies.

```
HP:0000455
HP:0000574
HP:0030084
HP:0012471
HP:0000239
HP:0001572
HP:0000960
HP:0001250
HP:0000322
HP:0001831
```

First copy one HPO term above to the search bar Phenomizer and click "Search".

![image](https://user-images.githubusercontent.com/11565618/173252411-dad15093-1ce0-4a19-942f-d3060c9714b5.png)

Then right click the searched HPO term and select "Add to Patient's Features". You will see this HPO term added to the right panel of Phenomizer. Repeat this step for all HPO terms listed above. 

![image](https://user-images.githubusercontent.com/11565618/173252457-7f675892-8a4b-4562-b617-2c10fbe3b351.png)

Select "Yes" to make sure 'symmetric' mode algorithm checked on. 

![image](https://user-images.githubusercontent.com/11565618/173253745-d951c9cd-6a64-4397-a98b-9306207a0aea.png)

After all HPO terms added to the right panel, click "Get diagnosis" at the bottom right corner of Phenomizer.

![image](https://user-images.githubusercontent.com/11565618/173253652-8b8714be-d676-43d5-8286-d7a5513c8f30.png)

You should see that KBG syndrome is among the top disease list

![image](https://user-images.githubusercontent.com/11565618/173252307-b87b7b5a-5de9-483a-a29c-d735a0760d09.png)

### 9. Running OARD
OARD (An acronym form of ["oh alright"](https://www.urbandictionary.com/define.php?term=oard)) is an **O**pen real-world based **A**nnotation for **R**are **D**iseases and its associated phenotypes.


This is a React web app to serve the web app of OARD. The backend is provided by OARD API. Currently it is hosted on the 
[NCATS AWS server (https://rare.cohd.io/)](https://rare.cohd.io/).

## An overview of OARD project
![image](https://github.com/WGLab/QuantitativeGenomics2023/assets/5926328/477c471c-a9e8-434b-b187-68993b122e3b)

### 1. Using OARD web app

Go to https://rare.cohd.io/.

![image](https://github.com/WGLab/QuantitativeGenomics2023/assets/5926328/54bb9097-2d4d-46a2-96b7-8004b09b9d07)

Select Method as `mostFrequency` and Domain as `Disease` to view most frequently occurred disease concept in a dataset:

![image](https://github.com/WGLab/QuantitativeGenomics2023/assets/5926328/8e6a4c0a-3c30-42ab-9403-1effedf34b98)

Then click Submit:

![image](https://github.com/WGLab/QuantitativeGenomics2023/assets/5926328/3585019f-7ec4-4622-9e70-7367f953cce1)

Type `HP:0012461` into Query Concept List 1 and select Method as `singleConceptFreq` and Dataset as `1-all` to view single concept frequency in a dataset:

![image](https://github.com/WGLab/QuantitativeGenomics2023/assets/11565618/0cd650f6-fba1-4221-b27b-89770574fade)

Then click Submit:

![image](https://github.com/WGLab/QuantitativeGenomics2023/assets/11565618/dce7beb6-eea3-4ec4-9c1b-de1344e7ee39)

Type `HP:0012461` into Query Concept List 1 and `HP:0500111` into Query Concept List 2 and select Method as `pairConceptFreq` and Dataset as `1-all` to view pair concept frequency in a dataset:

![image](https://github.com/WGLab/QuantitativeGenomics2023/assets/11565618/71f74440-013c-492a-b33c-f95baceca9f2)


### 2. Using OARD API
1. Go to Terminal, make sure you are in the `exercise2` directory first, and run 

`curl "https://rare.cohd.io/api/frequencies/mostFrequency?dataset_id=2&domain_id=diseases" > output_oard1.txt`

`curl "https://rare.cohd.io/api/frequencies/singleConceptFreq?dataset_id=1&concept_id=90012461" > output_oard2.txt`

`curl "https://rare.cohd.io/api/frequencies/pairedConceptFreq?dataset_id=1&concept_id_1=90012461&concept_id_2=90500111" > output_oard3.txt` 

where you generate JSON outputs.
However, since the `output.txt` file is in JSON format, it is not very intuitive to view the content of the file. Instead, we will use the table browser in Rstudio to view a formatted version of the JSON file.

2. Go To Console, remember that we are probably in the `exercise1` directory, so we should first set `exercise2` as the working directory.
```
setwd("../exercise2")
```
3. Then we can run
```
library("rjson")
# Read JSON results
result1 <- fromJSON(file = "output_oard1.txt")
# Convert them to array and name columns.
marray1 <- t(array( unlist(result1$results), dim=c(7, length(result1$results)) ) );
colnames(marray1) <- names(result1$results[[1]]);
# View the results in 2-D array. The second column is the rank of genes.
View (marray1);

result2 <- fromJSON(file = "output_oard2.txt")
# Convert them to array and name columns.
marray2 <- t(array( unlist(result2$results), dim=c(7, length(result2$results)) ) );
colnames(marray2) <- names(result2$results[[1]]);
# View the results in 2-D array. The second column is the rank of genes.
View (marray2);
```
![image](https://github.com/WGLab/QuantitativeGenomics2023/assets/5926328/a27c7a6d-c436-4384-8bbd-13b42df189d6)

### 10. Summary of the phenotype analysis exercises

In summary, a number of computational tools such as Phen2Gene, AMELIE and GADO can perform phenotype-driven gene prioritization. Phen2Gene provides webserver or API, or you can install and run locally (which is important to deploy it in batch processing mode), and it does not require a list of random genes to run either.

# CITATIONS:

- [ANNOVAR](https://academic.oup.com/nar/article/38/16/e164/1749458): Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data. Nucleic Acids Research, 38:e164, 2010
- [Phen2Gene](https://academic.oup.com/nargab/article/2/2/lqaa032/5843800): Zhao, M. et al. Phen2Gene: rapid phenotype-driven gene prioritization for rare diseases. NAR Genom Bioinform, 2:lqaa032 (2020).
- [AMELIE](https://stm.sciencemag.org/content/12/544/eaau9113.abstract): Birgmeier, J. et al. AMELIE speeds Mendelian diagnosis by matching patient phenotype and genotype to primary literature. Sci. Transl. Med. 12:eaau9113, (2020).
- [GADO](https://www.nature.com/articles/s41467-019-10649-4): Deelen, P. et al. Improving the diagnostic yield of exome- sequencing by predicting gene-phenotype associations using large-scale gene expression analysis. Nat. Commun. 10:2837 (2019).
- [ClinPhen](https://www.nature.com/articles/s41436-018-0381-1): Deisseroth, C. A. et al. ClinPhen extracts and prioritizes patient phenotypes directly from medical records to expedite genetic disease diagnosis. Genet. Med. 21:1585–1593 (2019).
- [PhenCards](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00909-8): Havrilla, J.M. et al. PhenCards: a data resource linking human phenotype information to biomedical knowledge. Genom. Med. 13:91 (2021)
- [OARD](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9502051/): Liu, C. et al. OARD: Open annotations for rare diseases and their phenotypes based on real-world data. Am. J. Hum. Genet. 109(9): 1591–1604, (2022)
