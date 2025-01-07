//----------------------------------------------------------------------------//
// TEXT TO DATA V3
// Purpose: Routine to take individual text files and places them in one free format text file that can be read as tabular data 
// Author : Steffen Erickson, Adapted from The Stata Journal (2009) 9, Number 4, pp. 599–620
// Date   : August 21, 2024 
//----------------------------------------------------------------------------//

mata 
struct myproblem {
	struct file_record scalar fr
	string scalar line
	real scalar output_fh
}
struct file_record {
	string scalar id 
	string scalar text
}
void driver(string scalar filespec, string scalar output_filename)
{
	string colvector filenames
	real scalar i
	real scalar output_fh
	filenames = sort(dir(".", "files", filespec),1)
	output_fh = fopen(output_filename, "w")
	for (i=1; i<=length(filenames); i++) {
		process_file(filenames[i], output_fh)
	}
	fclose(output_fh)
}
void process_file(string scalar filename, real scalar output_fh)
{
	struct myproblem scalar pr
	initialize_record(pr.fr)
	pr.output_fh = output_fh
	pr.fr.id = filename
	storetext(pr)
	output_record(pr)
}
void initialize_record(struct file_record scalar fr)
{
	fr.id = ""
	fr.text = ""
}
void storetext(struct myproblem scalar pr) 
{
	real scalar fh
	
	pr.fr.text = ""
	fh = fopen(pr.fr.id, "r")
		while ((pr.line=fget(fh))!=J(0,0,"")) {
				pr.fr.text = pr.fr.text + "\n" + pr.line 
			}
	fclose(fh)
	pr.fr.text = subinstr(pr.fr.text, `"""', "")
}
void output_record(struct myproblem scalar pr)
{
	fput(pr.output_fh, sprintf(`""%s" "%s""', pr.fr.id, pr.fr.text))
}
end 


