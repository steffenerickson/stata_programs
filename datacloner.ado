// Source https://www.stata.com/python/api16/Frame.html
python 
class DataCloner:
    def __init__(self, df, f):
	
	    #param df: The pandas DataFrame to clone data from.
        #param f: The target stata frame to clone variables and data into.

        self.df = df
        self.f = f

    def clone_variables(self):
        # Clone variables from the DataFrame to the target frame f 
		
        nvar = len(self.df.columns)

        for i in range(nvar):
            varname = Data.getVarName(i)
            vartype = Data.getVarType(i)
            if vartype == "byte":
                self.f.addVarByte(varname)
            elif vartype == "double":
                self.f.addVarDouble(varname)
            elif vartype == "float":
                self.f.addVarFloat(varname)
            elif vartype == "int":
                self.f.addVarInt(varname)
            elif vartype == "long":
                self.f.addVarLong(varname)
            elif vartype == "strL":
                self.f.addVarStrL(varname)
            else:
                self.f.addVarStr(varname, 10)

            self.f.setVarFormat(i, Data.getVarFormat(i))
            self.f.setVarLabel(i, Data.getVarLabel(i))

    def clone_data(self):

        # Clone data values from the DataFrame to the target frame f 
		
        self.f.setObsTotal(len(self.df))
        nvar = len(self.df.columns)

        for i in range(nvar):
            self.f.store(i, None, Data.get(var=i))
end 
