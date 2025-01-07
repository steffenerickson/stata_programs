// Source https://www.stata.com/python/api16/Frame.html



python 
def clone_data_to_stata(df, stata_frame):
    """
    Clone variables and data from a Pandas DataFrame to a Stata frame.
    """
    # Sanitize column names
    df.columns = (
        df.columns.str.strip()
        .str.replace(r"[^\w]", "_", regex=True)
        .str.replace(r"^\d", "var_", regex=True)
        .str[:32]
    )

    # Debug: Print sanitized column names
    print("Sanitized column names:", df.columns.tolist())

    # Clone variables
    for col in df.columns:
        vartype = df[col].dtype
        if pd.api.types.is_integer_dtype(vartype):
            stata_frame.addVarLong(col)
        elif pd.api.types.is_float_dtype(vartype):
            stata_frame.addVarDouble(col)
        elif pd.api.types.is_string_dtype(vartype):
            max_length = df[col].str.len().max() if not df[col].isnull().all() else 10
            stata_frame.addVarStr(col, max_length)
        else:
            raise ValueError(f"Unsupported data type for column: {col}")

        stata_frame.setVarLabel(stata_frame.getVarIndex(col), f"Label for {col}")

    # Clone data
    stata_frame.setObsTotal(len(df))
    for col in df.columns:
        col_index = stata_frame.getVarIndex(col)
        stata_frame.store(col_index, None, df[col].tolist())

end 



// Split into sentences
	* Requires pandas and nltk's sent_tokenize
	
python 
def split_by_sentence(df, text_column):
    # Check if index is single or multi-index
    if isinstance(df.index, pd.MultiIndex):
        index_names = df.index.names  # Use multi-index names
    else:
        index_names = [df.index.name or "index"]  # Use single-index name or default to "index"

    # Collect expanded rows
    expanded_rows = []
    for index_values, row in df.iterrows():
        # Ensure index_values is always a tuple for uniform processing
        if not isinstance(index_values, tuple):
            index_values = (index_values,)
        
        sentences = sent_tokenize(row[text_column])  # Split text into sentences
        for sentnum, sentence in enumerate(sentences, start=1):
            expanded_rows.append((*index_values, sentnum, sentence))  # Add sentence number

    # Create new DataFrame with expanded rows
    new_index_names = list(index_names) + ['sentnum']
    df_expanded = pd.DataFrame(
        expanded_rows, columns=new_index_names + [text_column]
    ).set_index(new_index_names)

    return df_expanded
end 


// tokenize sentences 
	* Requires pandas and nltk's word_tokenize and numpy
	
	
python 
def tokenize_sentences(df, text_column, keep_whitespace=True):
    # Ensure nltk resources are downloaded
    nltk.download('punkt', quiet=True)
    nltk.download('averaged_perceptron_tagger', quiet=True)

    # Tokenizer selection
    tokenizer = nltk.word_tokenize if keep_whitespace else nltk.WhitespaceTokenizer().tokenize

    # Process each sentence into tokens
    expanded_rows = []
    for index_values, row in df.iterrows():
        sentence = row[text_column]
        tokens = tokenizer(sentence)
        pos_tags = nltk.pos_tag(tokens)  # Part-of-speech tagging

        for tokennum, (token, pos) in enumerate(pos_tags, start=1):  # Add tokennum
            expanded_rows.append((*index_values, tokennum, token, pos))  # Retain index and add token info

    # Create a DataFrame for tokens
    token_index_names = list(df.index.names) + ['tokennum']
    df_tokens = pd.DataFrame(expanded_rows, columns=token_index_names + ['token_str', 'pos'])

    # Add additional columns for analysis
    df_tokens['term_str'] = df_tokens['token_str'].str.lower().str.replace(r"\W+", "", regex=True)
    df_tokens['pos_group'] = df_tokens['pos'].str[:2]

    return df_tokens
end 	
	
// Create BAG of words representation of the text
python 
def create_bow(CORPUS, bag, item_type='term_str'):
    BOW = CORPUS.groupby(bag+[item_type])[item_type].count().to_frame('n')
    return BOW
end 

// Calculate Term Frequency Inverse Document Frequency (TFIDF)
python 
from numpy.linalg import norm
from scipy.spatial.distance import pdist, squareform

def get_tfidf(BOW, tf_method='max', df_method='standard', item_type='term_str'):
            
    DTCM = BOW.n.unstack() # Create Doc-Term Count Matrix
    
    if tf_method == 'sum':
        TF = (DTCM.T / DTCM.T.sum()).T
    elif tf_method == 'max':
        TF = (DTCM.T / DTCM.T.max()).T
    elif tf_method == 'log':
        TF = (np.log2(DTCM.T + 1)).T
    elif tf_method == 'raw':
        TF = DTCM
    elif tf_method == 'bool':
        TF = DTCM.astype('bool').astype('int')
    else:
        raise ValueError(f"TF method {tf_method} not found.")

    DF = DTCM.count() # Assumes NULLs 
    N_docs = len(DTCM)
    
    if df_method == 'standard':
        IDF = np.log2(N_docs/DF) # This what the students were asked to use
    elif df_method == 'textbook':
        IDF = np.log2(N_docs/(DF + 1))
    elif df_method == 'sklearn':
        IDF = np.log2(N_docs/DF) + 1
    elif df_method == 'sklearn_smooth':
        IDF = np.log2((N_docs + 1)/(DF + 1)) + 1
    else:
        raise ValueError(f"DF method {df_method} not found.")
    
    TFIDF = TF * IDF
    
    DFIDF = DF * IDF
    
    TFIDF = TFIDF.fillna(0)

    return TFIDF, DFIDF
end 











