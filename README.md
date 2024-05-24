# deMix
## Hydrogen Deuterium eXchange Analysis 
Amide hydrogen/deuterium exchange (HDX) coupled with mass spectrometry (MS) is one of the most favorable tools for characterizing the protein dynamics and changes of protein conformation.
deMix is a fully automated software for the HDX-MS data analysis. deMix deals directly with the deuterated isotopic distributions, but not considering their centroid masses and is designed to be robust over random noises. 
In addition, deMix can also detect a bimodal deuterated distribution, arising from EX1 behavior or heterogeneous peptides in conformational isomer proteins.
<hr>


## 1. Usage
<pre>
- Command  : java -jar deMix.jar [options] [file]
  Options  :
  -i [parameter_file] : Config file path for search parameters [Required]
  -o [results_name]   : Name of analysis results [Optional]
                      : If not specified, the file name in 'Peptide=' parameter is used as the results_name.
                      : results_name_DdeuAnal.tsv and results_name_HDXProfile.tsv are generated.
- Example1 : java -jar deMix.jar -i foo.params -o hdx.foo
- Example2 : java -jar deMix.jar -i foo.params
- Requirement : deMix requires Java 1.7 or later. 
	        Visit http://www.oracle.com/technetwork/java/index.html 
</pre>
## 2. Parameters (See foo.txt file)
<pre>
- Peptide=[FILENAME]: 
  Specifies path to the file of peptides (in TSV format)
  The TSV file must include three columns for peptides, charge states, mz values.  
  Make sure that their column names (in the header) are 'peptide', 'charge', 'mz'. 
  They can be in any order. (See a sample file, foo_pept.tsv in the testdata directory, from https://prix.hanyang.ac.kr/download/deMixTestSample.zip)
  The TSV file may include different columns (they will be ignored).

- CTRLData= [FILENAME]: 
  Specifies the spectra file of a control experiment without HDX
  mzXML (centroid) file is allowed. 
	Elution time spans of input peptides are determined using this file.

- HDXData=[Label],[spectra file]: 
  Specifies the spectra file of HDX experiments with the D2O labeling
  mzXML (centroid) file is allowed. 
  Multiple files can be specified using different labeling as below.
  e.g., HDXData=30s, D:/mydir/my_d2o_30s.mzXML
        HDXData=20m, D:/mydir/my_d2o_10m.mzXML
        HDXData=30m, D:/mydir/my_d2o_30m.mzXML

- MassTolerance=[MASS_TOLERANCE]: 
  Sets a peptide mass tolerance in dalton or ppm. The default value is 10ppm.
  e.g., MassTolerance=15ppm
  e.g., MassTolerance=0.4da

- Protein=[FILENAME]: 
  Specifies path to the file of protein sequences (*.fasta format)
  Optional. If it's specified, start and end positons of peptides will be reported.
  e.g., Protein= D:/mydir/my_protein.fasta
</pre>
## 3. Citation
<pre>
- deMix: Decoding deuterated distributions from heterogeneous protein states via HDX-MS.
  Na S, Lee JJ, Joo JWJ, Lee KJ, Paek E. 
  Scientific Reports <b>2019</b>, 9(1), 3176.
</pre>
## 4. Rights and Permissions
<pre>
- deMix © 2024 is licensed under Creative Commons Attribution-ShareAlike 4.0 International.
  This license requires that reusers give credit to the creator. It allows reusers to distribute, 
  remix, adapt, and build upon the material in any medium or format, even for commercial purposes. 
  If others remix, adapt, or build upon the material, they must license the modified material under identical terms.
  To view a copy of this license, visit https://creativecommons.org/licenses/by-sa/4.0/
</pre>

