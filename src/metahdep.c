
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <stdio.h>
#include <string.h>			//	strcpy, strcasecmp


/*******************************************************************************************

	metan_merge :  merge two sorted sub-arrays into one sorted subarray
				   this function is needed by the mergesort function

	               name_list is the SEXP with all the gene names
				   sort_index is an int pointer to hold the sorted indices
				   left_part[] and right_part[] are the subarrays holding the 2 index lists
				   llen and rlen and the lengths of the left and right parts

  ******************************************************************************************/
void metan_merge(SEXP name_list, int sort_index[], int left_part[], int llen, int right_part[], int rlen)
{
	int i=0,
		li=0,
		ri=0;

	while (li < llen && ri < rlen)
	{
		//	compare the gene names pointed to by the two indices and merge accordingly
		//	strcasecmp ignores case
		if (strcasecmp(CHAR(STRING_ELT(name_list, left_part[li])),
			       CHAR(STRING_ELT(name_list, right_part[ri]))) <=0)
			sort_index[i++] = left_part[li++];
		else
			sort_index[i++] = right_part[ri++];
	}

	//	append the remaining indices from whichever sub-array still has some left
	while (li < llen)
		sort_index[i++] = left_part[li++];
	while (ri < rlen)
		sort_index[i++] = right_part[ri++];
}


/*******************************************************************************************

	metan_mergesort :  recursively sort an array of indices which when used alphabetize the
					   name_list vector of strings

					   name_list is the list of gene names in R
					   sort_indices is an integer array holding the indices to the names
					   first and last are the indices first and last elements to be sorted
					   
 *******************************************************************************************/
void metan_mergesort(SEXP name_list, int sort_indices[], int first, int last)
{
	int i,
		middle,
		*temp_array;

	//	if the section of the array hasn't been fully split,
	//	then split it and (recursively) sort each side
	if (first != last)
	{
		//	split the array in half and sort each half
		middle = first + (last - first) / 2;
		metan_mergesort(name_list, sort_indices, first, middle);
		metan_mergesort(name_list, sort_indices, middle+1, last);
	
		//	then, merge the two sides into a new array
		temp_array = (int *) R_alloc((last - first + 1), sizeof(int));
		metan_merge(name_list, temp_array, &(sort_indices[first]), middle-first+1,
										   &(sort_indices[middle+1]), last-middle);
		
		//	copy the contents of the new array into place in the original array
		for (i=0; i<last-first+1; i++)
			sort_indices[first+i] = temp_array[i];
	}
}


/*******************************************************************************************

	sort_newname_list :  create an index for a sorted newname list, sorted by new.name
						 unlike the old gene names, which all appear to be unique within
						 each study, the new.names are not unique and there will be multiple
						 rows having the same new.name

						 returns a vector of integers containing the sorted indices

  ******************************************************************************************/
SEXP sort_newname_list(SEXP gene_list)
{
	int i,					//	generic loop index
		*temp_ptr,			//	temporary ptr to the R vector to be returned
		num_names,			//	number of names in the passed vector gene_list
		*sort_index;		//	array of indices to be sorted

	SEXP return_vector;		//	R structure for returning the sorted index array


	//	specify that this is a vector of character strings
	PROTECT(gene_list = AS_CHARACTER(gene_list));
	num_names = length(gene_list);

	//	allocate memory for the sorted indices and initialize them
	sort_index = (int *) R_alloc(num_names, sizeof(int));
	for (i=0; i<num_names; i++)
		sort_index[i] = i;

	//	sort the indices using mergesort
	metan_mergesort(gene_list, sort_index, 0, num_names-1);

	//	create the return vector and copy the sorted indices into it
	PROTECT(return_vector = NEW_INTEGER(num_names));
	temp_ptr = INTEGER_POINTER(return_vector);
	//memcpy(temp_ptr, sort_index, num_names*sizeof(int));
	for (i=0; i<num_names; i++)
		temp_ptr[i] = sort_index[i];

	UNPROTECT(2);	//	gene_list, return_vector

	return(return_vector);
}


/*******************************************************************************************

	sort_name_lists :  create an index for each vector of passed gene names.  Sort the indices
					   so that by taking the gene names in the order given in the indices then
					   the genes will be ordered alphabetically.  i.e. the actual gene name
					   vectors are left unsorted; but the indices can be used to easily get the
					   ordered names.

					   returns a list of integer vectors

  ******************************************************************************************/
SEXP sort_name_lists(SEXP R_list_of_name_vectors)
{
	int	i, j,
		*temp_ptr,
		num_studies=0,
		study_length=0,
		**index_ptr_array;

	SEXP temp_vector,
		 return_list = R_NilValue,
		 *list_elements;

	//	need to get the number of studies and how many names are in each
	//	need to allocate temporary storage for the indices
	//	this is complicated by the fact that there can be a varying number of
	//	studies, each with a different number of genes
	num_studies = length(R_list_of_name_vectors);
	index_ptr_array = (int**) R_alloc(num_studies, sizeof(int *));
	for (i=0; i<num_studies; i++)
	{
		//	extract the 'i'th vector of names from the R list object
		PROTECT(temp_vector = AS_CHARACTER(VECTOR_ELT(R_list_of_name_vectors, i)));
		study_length = length(temp_vector);

		//	allocate space for the indices and initialize them
		index_ptr_array[i] = (int *) R_alloc(study_length, sizeof(int));
		for (j=0; j<study_length; j++)
		{  index_ptr_array[i][j] = j;  }

		//	perform a mergesort on the indices for the current study
		metan_mergesort(temp_vector, index_ptr_array[i], 0, study_length-1);

		UNPROTECT(1);	//	temp_vector
	}

	//	create R objects to hold the vectors of sorted indices, then copy the data into them
	list_elements = (SEXP *) R_alloc(num_studies, sizeof(SEXP));
	for (i=0; i<num_studies; i++)
	{
		study_length = length(VECTOR_ELT(R_list_of_name_vectors, i));
		PROTECT(list_elements[i] = NEW_INTEGER(study_length));
		temp_ptr = INTEGER_POINTER(list_elements[i]);
		for (j=0; j<study_length; j++)
			temp_ptr[j] = index_ptr_array[i][j];
	}
	
	//	alloc space for the return_list and insert the vectors
	PROTECT(return_list = allocVector(VECSXP, num_studies));
	for (i=0; i<num_studies; i++)
		SET_VECTOR_ELT(return_list, i, list_elements[i]);
	UNPROTECT(num_studies);	//	list_elements[]
	
	UNPROTECT(1);	//	return_list
	return(return_list);
}


/*******************************************************************************************

	metan_binary_search_not_unique :  performs a binary search via the passed array of
				sorted indices.  The not_unique refers to the fact that the gene name being
				searched for may appear multiple times in the name_list (but since we can
				treat the list as being sorted they will all be grouped together.)

				returns a length-2 vector of integers, being the first and last indices in
				the sort_index of the group of indices corresponding to gene_name

  ******************************************************************************************/
SEXP metan_binary_search_not_unique(SEXP name_list, SEXP sort_index, SEXP gene_name)
{
	int found = 0;

	int	first,
		middle=0,
		last,
		num_names,
		comp_result,
		*sort_index_ptr,
		*return_vector_ptr;

	char *cur_name;

	SEXP return_vector;

	//	extract the gene name into a more C-friendly form
	cur_name = (char *) R_alloc(strlen(CHAR(STRING_ELT(gene_name, 0)))+1, sizeof(char));
	strcpy(cur_name, CHAR(STRING_ELT(gene_name, 0)));

	//	get an int pointer to the sort index for easy C access
	sort_index_ptr = INTEGER_POINTER(sort_index);
	
	num_names = length(name_list);
	first = 0;
	last = num_names - 1;

	//	do a binary search using the sort indices
	while (!found && first<=last)
	{
		middle = first + (last - first) / 2;
		comp_result = strcasecmp(cur_name, CHAR(STRING_ELT(name_list, sort_index_ptr[middle])));
		if (!comp_result)
			found = 1;
		else
		{
			if (comp_result < 0)
				last = middle - 1;
			else
				first = middle + 1;
		}
	}

	//	if the gene name was not found, return the R_NilValue so the calling function
	//	can respond appropriately
	if (!found)
		return(R_NilValue);

	//	at this point, all we know is that the gene name pointed to by 'middle' matches the cur_name
	//	we need to find the full range of indices that match cur_name
	first = middle;
	last = middle;
	
	//	if the gene in the spot before the one we found also matches, then include that one also
	//	keep checking each previous name in the list until a non-matching gene name is found
	//	this may be a little more elaborate than it needs to be
	//	it is done this way to prevent indexing into out-of-bounds areas
	//	i.e. before the beginning or after the end of the name list
	while (TRUE)
	{
		if (first - 1 < 0)
			break;
		else
		{
			if (!strcasecmp(cur_name, CHAR(STRING_ELT(name_list, sort_index_ptr[first - 1]))))
				first--;
			else
				break;
		}
	}

	//	do the same thing only with the subsequent gene instead of the previous gene
	while (TRUE)
	{
		if (last + 1 > num_names - 1)
			break;
		else
		{
			if (!strcasecmp(cur_name, CHAR(STRING_ELT(name_list, sort_index_ptr[last + 1]))))
				last++;
			else
				break;
		}
	}

	//	create a new object to hold the data to return
	PROTECT(return_vector = NEW_INTEGER(2));
	return_vector_ptr = INTEGER_POINTER(return_vector);
	return_vector_ptr[0] = first;
	return_vector_ptr[1] = last;

	UNPROTECT(1);	//	return_vector
	return(return_vector);
}


/*******************************************************************************************

	find_old_names2 :  looks through the name_list for all of the indices corresponding to
			records having the same new.name (Unigene ID, etc.) as the passed gene_name

			returns a vector of integers that can be used to create a list of old.names
			(probeset names) and chipset names, which must then be found in the main data
			matrix.

			the indices refer to rows in the newnames dataframe.  the old.name on the
			referenced rows are the old.names corresponding the the passed new.name

  ******************************************************************************************/
SEXP find_old_names2(SEXP name_list, SEXP sort_index, SEXP gene_name)
{
	int i,						//	loop counter
		first,					//	the first index of a section of the array
		last,					//	the last index of a section of the array
		*sort_index_ptr,		//	a pointer used to copy the data into the return_list
		*return_list_ptr;		//	another pointer used to copy data into the return_list

	SEXP range_indices,
		 return_list;

	//	do a binary search using the sort indices
	PROTECT(range_indices = metan_binary_search_not_unique(name_list, sort_index, gene_name));

	//	if the gene is not found, then return R_NilValue so the calling function
	//	(this will be called via .Call from R) can respond appropriately
	if (range_indices == R_NilValue)
	{
		UNPROTECT(1);	//	range_indices
		return(R_NilValue);
	}

	//	get the range of sort_index indices corresponding to the gene in gene_name
	first = INTEGER_POINTER(range_indices)[0];
	last = INTEGER_POINTER(range_indices)[1];
	UNPROTECT(1);	//	range_indices

	//	expand the index list and insert it into an R vector
	PROTECT(return_list = NEW_INTEGER(last - first + 1));
	return_list_ptr = INTEGER_POINTER(return_list);
	sort_index_ptr = INTEGER_POINTER(sort_index);

	return_list_ptr[0] = sort_index_ptr[first];
	for (i=first; i<=last; i++)
		return_list_ptr[i-first] = sort_index_ptr[i];

	UNPROTECT(1);	//	return_list
	return(return_list);
}


/*******************************************************************************************

	metan_binary_search_unique :  performs a binary search via the passed array of sorted indices.

				the _unique indicates that the function only expects one find one instance of a
				gene name in the passed name list, i.e. there shouldn't be multiple rows with the
				same gene name.  if there are, it still only returns the row index of the first
				one it runs across, so it may be unpredictable if the data frames are not set up
				correctly.  the function will still return a valid index(or R_NilValue) but *which*
				row index (in the case of multiple rows with the same name) will change depending on
				the number of genes in the list.

				returns the index of the row in which the match occurred, or R_NilValue if
				the gene name was not found in the list

  ******************************************************************************************/
SEXP metan_binary_search_unique(SEXP name_list, SEXP sort_index, SEXP gene_name, int gene_index, char chipset_name, int study_num)
{
	int found = 0;

	int	first,
		middle=0,
		last,
		num_names,
		comp_result,
		*sort_index_ptr,
		*return_index_ptr;

	char *cur_name;

	SEXP return_index;


	//	extract the gene name into a more C-friendly form
	PROTECT(gene_name = AS_CHARACTER(gene_name));
	cur_name = (char *) R_alloc(strlen(CHAR(STRING_ELT(gene_name, gene_index))), sizeof(char));
	strcpy(cur_name, CHAR(STRING_ELT(gene_name, gene_index)));
	UNPROTECT(1);	//	gene_name

	//	get an int pointer to the sort index for easy C access
	sort_index_ptr = INTEGER_POINTER(sort_index);
	num_names = length(name_list);
	first = 0;
	last = num_names - 1;
	
	while (!found && first<=last)
	{
		middle = first + (last - first) / 2;
		comp_result = strcasecmp(cur_name, CHAR(STRING_ELT(name_list, sort_index_ptr[middle])));

		if (!comp_result)
			found = 1;
		else
		{
			if (comp_result < 0)
				last = middle - 1;
			else
				first = middle + 1;
		}
	}

	//	if the gene name was not found, then return R_NilValue so the
	//	calling function can respond appropriately
	if (!found)
		return(R_NilValue);

	PROTECT(return_index = NEW_INTEGER(1));
	return_index_ptr = INTEGER_POINTER(return_index);
	return_index_ptr[0] = middle;

	UNPROTECT(1);	//	return_index
	return(return_index);
}


#define MAX_NUM_ROWS 300
SEXP get_row_indices2(SEXP R_study_gene_names, SEXP R_old_names, SEXP R_old_chipsets, SEXP R_study_chipsets, SEXP full_sort_index)
{
	int	i, j, k,				//	loop counters
		num_rows_found = 0,		//	the total number of records found for the current gene (accumulates)
		num_names,				//	the number of old.names that need to be matched in the study arrays
		num_studies,			//	the number of studies in the passed data list
		*study,					//	a parallel array to hold the study numbers of found genes
		*row,					//	a parallel array to hold the row numbers of found genes
		*temp1, *temp2,			//	used when expanding the number of allowable rows
		*ireturn_matrix,		//	a pointer used to fill the R return matrix
	    max_rows=MAX_NUM_ROWS;	//	the current maximum number of rows

	SEXP index,					//	object to hold the results of the search
		 name_list,
		 sort_index,
		 return_matrix;			//	object to hold the matrix of returned study-row paired indices


	//	allocate space to store the row information
	study = (int *) R_alloc(max_rows, sizeof(int));
	row = (int *) R_alloc(max_rows, sizeof(int));

	//	extract numbers to avoid function calls in the middle of loops (just in case)
	num_studies = length(R_study_gene_names);
	num_names = length(R_old_names);

	//	loop through the different studies
	//	for each study, find get the chipset-specific gene/probeset names
	//	corresponding to the current gene of interest
	for (i=0; i<num_studies; i++)
	{
		for (j=0; j<num_names; j++)
		{
			if (!strcasecmp(CHAR(STRING_ELT(R_study_chipsets, i)), CHAR(STRING_ELT(R_old_chipsets, j))))
			{
				name_list = VECTOR_ELT(R_study_gene_names, i);
				sort_index = VECTOR_ELT(full_sort_index, i);

				//	do a binary search for a SINGLE row
				index = metan_binary_search_unique(name_list, sort_index, R_old_names, j, CHAR(STRING_ELT(R_old_chipsets, j)), i);

				//	make sure the name was found before adding it to the row.indices
				if (index != R_NilValue)
				{
					//	check to see if there is room for the new row
					//	if not, allocate more space and copy the contents into the larger spot
					if (num_rows_found + 1 > max_rows)
					{
						max_rows += MAX_NUM_ROWS;
						temp1 = study;
						temp2 = row;
						study = (int *) R_alloc(max_rows, sizeof(int));
						row = (int *) R_alloc(max_rows, sizeof(int));
						for (k=0; k<num_rows_found; k++)
						{
							study[k] = temp1[k];
							row[k] = temp2[k];
						}
					}

					//	extract the index and save the study and row numbers
					PROTECT(index = AS_INTEGER(index));
					study[num_rows_found] = i;
					row[num_rows_found] = *INTEGER_POINTER(index);
					num_rows_found++;
					UNPROTECT(1);	//	index
				}
			}
		}
	}

	//	if no rows were found then return R_NilValue so the calling function knows
	if (!num_rows_found)
		return(R_NilValue);

	//	insert the results into a nice matrix and return it
	PROTECT(return_matrix = allocMatrix(INTSXP, num_rows_found, 2));
	ireturn_matrix = INTEGER(return_matrix);
	for (i=0; i<num_rows_found; i++)
	{
		ireturn_matrix[i] = study[i];
		ireturn_matrix[num_rows_found + i] = row[i];
	}

	UNPROTECT(1);	//	return_matrix
	return(return_matrix);
}


