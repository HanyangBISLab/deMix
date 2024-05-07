# deMix
## Hydrogen Deuterium eXchange Analysis 

<hr>


## 1. Usage
<pre>
- Command  : java -jar deMix.jar [options] [file]
  Options  :
  -i <parameter_file> : Config file path for search parameters [Required]
  -o <results_name>   : Name of analysis results [Optional]
                      : If not specified, the file name in 'Peptide=' parameter is used as the results_name.
                      : results_name_DdeuAnal.tsv and results_name_HDXProfile.tsv are generated.
- Example1 : java -jar deMix.jar -i foo.params -o hdx.foo
- Example2 : java -jar deMix.jar -i foo.params
- Requirement : deMix requires Java 1.7 or later. Visit http://www.oracle.com/technetwork/java/index.html 
</pre>
## 2. Parameters (See foo.txt file)
<pre>

</pre>
## 3. Citation
<pre>
- deMix: Decoding deuterated distributions from heterogeneous protein states via HDX-MS.
  Na S, Lee JJ, Joo JWJ, Lee KJ, Paek E. 
  Scientific Reports <b>2019</b>, 9(1), 3176.
</pre>
## 4. Rights and Permissions
<pre>
- deMix by 2024 is licensed under Creative Commons Attribution-ShareAlike 4.0 International, 
  which permits use, sharing, adaptation, distribution and reproduction in any medium or format, 
  as long as you give appropriate credit to the original author(s) and the source, 
  provide a link to the Creative Commons license, and indicate if changes were made.
  To view a copy of this license, visit https://creativecommons.org/licenses/by-sa/4.0/
</pre>

