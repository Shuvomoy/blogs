@def title = "Working with dataframes in Julia"
@def published ="May 5, 2020"
@def tags =["programming", "Julia"]


# Working with dataframes in Julia
**Shuvomoy Das Gupta**

*May 5, 2020*




In this blog, we will discuss how to work with dataframes using the `DataFrames` package in `Julia`. As an illustrative example, we will work with the MovieLens dataset. This is an introductory blog, and for learning how to use the `DataFrames` package in greater details, a great set of tutorials is available with jupyter notebooks at [this link](https://github.com/bkamins/Julia-DataFrames-Tutorial/).

---

**Table of contents**

\toc

---

**Jupyter notebook for this blog.** The jupyter notebook for this blog can be downloaded from [this link](https://raw.githubusercontent.com/Shuvomoy/blog/gh-pages/codes/working_with_dataframes_intro_Julia.ipynb) and the data files in zip format are available [here](https://github.com/Shuvomoy/blog/blob/gh-pages/data/movielens_datasets.zip). Please unzip them in the folder, where your ipynb or julia file is.

**What is the MovieLens dataset?** The `MovieLens` dataset is one of the most common datasets people use for testing recommendation system algorithm. The version I will work with contains 1,000,209 ratings for approximately 3,900 movies; these recommendations were made by 6,040 MovieLens users. 

Let us start with the necessary packages. 


```julia
using CSV, DataFrames
```

First, let us check if the files are loaded correctly. The files are available here: **Provide the link**.


```julia
# First, let us check if the files are loaded correctly
cd("C:\\Users\\shuvo\\Desktop\\MovieLensDataSet_Experiment")
isfile.(["users.csv" "movies.csv" "ratings.csv"])
```




    1×3 BitArray{2}:
     1  1  1



### Putting the data into dataframes
Let us put the data into dataframes.


```julia
# put the data in dataframes
df_users = CSV.read("users.csv");
df_movies = CSV.read("movies.csv");
df_ratings = CSV.read("ratings.csv");
```

**The `df_users` dataframe.** Let us briefly explore the contents of the `df_users` dataframe. 


```julia
size(df_users)
```




    (6040, 8)



So, `df_users` has 6040 rows and 8 columns. It is too large, but let us peek into the first and last few of the rows rows.


```julia
first(df_users,3)
```

~~~


<table class="data-frame"><thead><tr><th></th><th>Column1</th><th>user_id</th><th>gender</th><th>age</th><th>occupation</th><th>zipcode</th><th>age_desc</th><th>occ_desc</th></tr><tr><th></th><th>Int64</th><th>Int64</th><th>String</th><th>Int64</th><th>Int64</th><th>String</th><th>String</th><th>String</th></tr></thead><tbody><p>3 rows × 8 columns</p><tr><th>1</th><td>0</td><td>1</td><td>F</td><td>1</td><td>10</td><td>48067</td><td>Under 18</td><td>K-12 student</td></tr><tr><th>2</th><td>1</td><td>2</td><td>M</td><td>56</td><td>16</td><td>70072</td><td>56+</td><td>self-employed</td></tr><tr><th>3</th><td>2</td><td>3</td><td>M</td><td>25</td><td>15</td><td>55117</td><td>25-34</td><td>scientist</td></tr></tbody></table>


~~~


```julia
last(df_users, 4)
```

~~~



<table class="data-frame"><thead><tr><th></th><th>Column1</th><th>user_id</th><th>gender</th><th>age</th><th>occupation</th><th>zipcode</th><th>age_desc</th><th>occ_desc</th></tr><tr><th></th><th>Int64</th><th>Int64</th><th>String</th><th>Int64</th><th>Int64</th><th>String</th><th>String</th><th>String</th></tr></thead><tbody><p>4 rows × 8 columns</p><tr><th>1</th><td>6036</td><td>6037</td><td>F</td><td>45</td><td>1</td><td>76006</td><td>45-49</td><td>academic/educator</td></tr><tr><th>2</th><td>6037</td><td>6038</td><td>F</td><td>56</td><td>1</td><td>14706</td><td>56+</td><td>academic/educator</td></tr><tr><th>3</th><td>6038</td><td>6039</td><td>F</td><td>45</td><td>0</td><td>01060</td><td>45-49</td><td>other or not specified</td></tr><tr><th>4</th><td>6039</td><td>6040</td><td>M</td><td>25</td><td>6</td><td>11106</td><td>25-34</td><td>doctor/health care</td></tr></tbody></table>

~~~

Now let us take a look at the basic summary statistics of data in df_users.


```julia
# Let us take a look at a quick summary about df_users for a few users
describe(df_users, cols = 1:3)
# if we were interested about all the columns, then we could run the following command:
# describe(df_users)
```



~~~
<table class="data-frame"><thead><tr><th></th><th>variable</th><th>mean</th><th>min</th><th>median</th><th>max</th><th>nunique</th><th>nmissing</th><th>eltype</th></tr><tr><th></th><th>Symbol</th><th>Union…</th><th>Any</th><th>Union…</th><th>Any</th><th>Union…</th><th>Nothing</th><th>DataType</th></tr></thead><tbody><p>3 rows × 8 columns</p><tr><th>1</th><td>Column1</td><td>3019.5</td><td>0</td><td>3019.5</td><td>6039</td><td></td><td></td><td>Int64</td></tr><tr><th>2</th><td>user_id</td><td>3020.5</td><td>1</td><td>3020.5</td><td>6040</td><td></td><td></td><td>Int64</td></tr><tr><th>3</th><td>gender</td><td></td><td>F</td><td></td><td>M</td><td>2</td><td></td><td>String</td></tr></tbody></table>
~~~


The get the name of all the columns of `df_users`, we can use the function `names`. 


```julia
names(df_users)
```




    8-element Array{Symbol,1}:
     :Column1   
     :user_id   
     :gender    
     :age       
     :occupation
     :zipcode   
     :age_desc  
     :occ_desc  



If we want to know the element types of each column, we can run the following.


```julia
eltype.(eachcol(df_users))
```




    8-element Array{DataType,1}:
     Int64 
     Int64 
     String
     Int64 
     Int64 
     String
     String
     String



**The `df_movies` dataframe.** Let us repeat the previous commands for the `df_movies` dataframe.


```julia
size(df_movies);
```


```julia
first(df_movies,5);
```


```julia
last(df_movies, 3);
```


```julia
describe(df_movies, cols = 1:3);
```


```julia
names(df_movies)
```




    4-element Array{Symbol,1}:
     :Column1 
     :movie_id
     :title   
     :genres  




```julia
eltype.(eachcol(df_movies));
```

**The `df_ratings` dataframe.** Let us see what is going on with the dataframe `df_ratings`.


```julia
size(df_ratings); # quite large with 1000209 columns and 7 rows
```


```julia
first(df_ratings,3)
```



~~~
<table class="data-frame"><thead><tr><th></th><th>Column1</th><th>user_id</th><th>movie_id</th><th>rating</th><th>timestamp</th><th>user_emb_id</th><th>movie_emb_id</th></tr><tr><th></th><th>Int64</th><th>Int64</th><th>Int64</th><th>Int64</th><th>Int64</th><th>Int64</th><th>Int64</th></tr></thead><tbody><p>3 rows × 7 columns</p><tr><th>1</th><td>0</td><td>1</td><td>1193</td><td>5</td><td>978300760</td><td>0</td><td>1192</td></tr><tr><th>2</th><td>1</td><td>1</td><td>661</td><td>3</td><td>978302109</td><td>0</td><td>660</td></tr><tr><th>3</th><td>2</td><td>1</td><td>914</td><td>3</td><td>978301968</td><td>0</td><td>913</td></tr></tbody></table>
~~~



```julia
last(df_ratings,3);
```


```julia
describe(df_ratings, cols = 1:3);
```


```julia
names(df_ratings);
```


```julia
eltype.(eachcol(df_ratings));
```


```julia
## Let us create the rating matrix
# first let us create the number of unique users and matrices
n_users = length(unique(df_ratings.user_id))
n_movies = length(unique(df_ratings.movie_id))
println("Number of users = $(n_users) | Number of movies = $(n_movies)")
```

    Number of users = 6040 | Number of movies = 3706


### Converting a dataframe into wide format
Now we construct a wide format dataframe for the `df_ratings` dataframe using the `unstack` function. Our goal is creating a dataframe, which will have the users as the rows, the movies as the columns, and the rating as the value. The `unstack` function requires specifying which columns would act as id variable, column variable name, and column values. 


```julia
df_ratings_unstacked = unstack(df_ratings, :user_id, :movie_id, :rating);
```


```julia
size(df_ratings_unstacked)
```




    (6040, 3707)




```julia
first(df_ratings_unstacked, 6)
```



~~~
<table class="data-frame"><thead><tr><th></th><th>user_id</th><th>1</th><th>2</th><th>3</th><th>4</th><th>5</th><th>6</th><th>7</th><th>8</th></tr><tr><th></th><th>Int64</th><th>Int64⍰</th><th>Int64⍰</th><th>Int64⍰</th><th>Int64⍰</th><th>Int64⍰</th><th>Int64⍰</th><th>Int64⍰</th><th>Int64⍰</th></tr></thead><tbody><p>6 rows × 3,707 columns (omitted printing of 3698 columns)</p><tr><th>1</th><td>1</td><td>5</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td></tr><tr><th>2</th><td>2</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td></tr><tr><th>3</th><td>3</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td></tr><tr><th>4</th><td>4</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td></tr><tr><th>5</th><td>5</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>2</td><td>missing</td><td>missing</td></tr><tr><th>6</th><td>6</td><td>4</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td><td>missing</td></tr></tbody></table>
~~~


### Working with missing values
We have lot of missing values, which we can see the let us replace the missing values with `NaN`.


```julia
# replace the missing values in the data frame
column_names_ru = names(df_ratings_unstacked);
for i in 1:length(column_names_ru)
   df_ratings_unstacked[!, column_names_ru[i]] = replace(df_ratings_unstacked[!, column_names_ru[i]], missing => NaN) # If we want to replace the missing values with any other number, then we use: missing => NaN
end

```

Let us see if the values have indeed changed.


```julia
first(df_ratings_unstacked, 6)

```



~~~
<table class="data-frame"><thead><tr><th></th><th>user_id</th><th>1</th><th>2</th><th>3</th><th>4</th><th>5</th><th>6</th><th>7</th><th>8</th></tr><tr><th></th><th>Float64</th><th>Float64</th><th>Float64</th><th>Float64</th><th>Float64</th><th>Float64</th><th>Float64</th><th>Float64</th><th>Float64</th></tr></thead><tbody><p>6 rows × 3,707 columns (omitted printing of 3698 columns)</p><tr><th>1</th><td>1.0</td><td>5.0</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td></tr><tr><th>2</th><td>2.0</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td></tr><tr><th>3</th><td>3.0</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td></tr><tr><th>4</th><td>4.0</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td></tr><tr><th>5</th><td>5.0</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>2.0</td><td>NaN</td><td>NaN</td></tr><tr><th>6</th><td>6.0</td><td>4.0</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td></tr></tbody></table>
~~~


**What if the missing values are represented as some other symbol in the original dataframe?** In the previous example, the missing values were convenienetly recognized by Julia as `missing` type. Now, in some data files the missing values are represented by some other symbols, such as "?". In such case the column that contains "?" may be diagnosed as an array of strings by Julia. In that case, we could run the following command to replace the "?"s with `missing` using the following command.


```julia
type_of_columns = eltype.(eachcol(df)) # say the name of the dataframe is df

m, n = size(df)

for i in 3:n # let the important columns containing missing values start from column number 3

  if type_of_columns[i] == String # due to presence of "?" in the column, it may get diagnosed by an array of strings by Julia
  
          df[!, i]  = map(x->
                         begin
                         val = tryparse(Float64, x)
                         ifelse(typeof(val) == Nothing, missing, val)                         
                         end, df[!,i])
                         
  end
end
```

### Saving and loading dataframe
Let us see how to load and safe a dataframe. We can save `df_ratings_unstacked` using the `write` function provided by the `CSV` package. 


```julia
CSV.write("df_ratings_unstacked.csv", df_ratings_unstacked)
```




    "df_ratings_unstacked.csv"



To load back the dataframe that we just saved into the csv file, we can use `read` command provided by `CSV` package. The command `use_mmap = false` is used, so that the file can be deleted in the same session; this is required only in Windows.


```julia
df_loaded = CSV.read("df_ratings_unstacked.csv", use_mmap=false);
first(df_loaded, 5)
```



~~~
<table class="data-frame"><thead><tr><th></th><th>user_id</th><th>1</th><th>2</th><th>3</th><th>4</th><th>5</th><th>6</th><th>7</th><th>8</th></tr><tr><th></th><th>Float64</th><th>Float64</th><th>Float64</th><th>Float64</th><th>Float64</th><th>Float64</th><th>Float64</th><th>Float64</th><th>Float64</th></tr></thead><tbody><p>5 rows × 3,707 columns (omitted printing of 3698 columns)</p><tr><th>1</th><td>1.0</td><td>5.0</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td></tr><tr><th>2</th><td>2.0</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td></tr><tr><th>3</th><td>3.0</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td></tr><tr><th>4</th><td>4.0</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td></tr><tr><th>5</th><td>5.0</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>NaN</td><td>2.0</td><td>NaN</td><td>NaN</td></tr></tbody></table>
~~~


### Converting a dataframe into a matrix
Now let us convert the `df_ratings_unstacked` dataframe into a matrix.


```julia
ratings_matrix = df_ratings_unstacked[:,  2:length(column_names_ru)];
rating_matrix = convert(Matrix, ratings_matrix)

```



```julia 
    6040×3706 Array{Float64,2}:
       5.0  NaN    NaN    NaN    NaN    …  NaN  NaN  NaN    NaN  NaN  NaN  NaN
     NaN    NaN    NaN    NaN    NaN       NaN  NaN  NaN    NaN  NaN  NaN  NaN
     NaN    NaN    NaN    NaN    NaN       NaN  NaN  NaN    NaN  NaN  NaN  NaN
     NaN    NaN    NaN    NaN    NaN       NaN  NaN  NaN    NaN  NaN  NaN  NaN
     NaN    NaN    NaN    NaN    NaN       NaN  NaN  NaN    NaN  NaN  NaN  NaN
       4.0  NaN    NaN    NaN    NaN    …  NaN  NaN  NaN    NaN  NaN  NaN  NaN
     NaN    NaN    NaN    NaN    NaN       NaN  NaN  NaN    NaN  NaN  NaN  NaN
       4.0  NaN    NaN      3.0  NaN       NaN  NaN  NaN    NaN  NaN  NaN  NaN
       5.0  NaN    NaN    NaN    NaN       NaN  NaN    3.0  NaN  NaN  NaN  NaN
       5.0    5.0  NaN    NaN    NaN       NaN  NaN    4.0  NaN  NaN  NaN  NaN
     NaN    NaN    NaN    NaN    NaN    …  NaN  NaN  NaN    NaN  NaN  NaN  NaN
     NaN    NaN    NaN    NaN    NaN       NaN  NaN  NaN    NaN  NaN  NaN  NaN
     NaN      3.0  NaN    NaN    NaN       NaN  NaN  NaN    NaN  NaN  NaN  NaN
       ⋮                                ⋱         ⋮                          ⋮
     NaN    NaN    NaN    NaN    NaN       NaN  NaN  NaN    NaN  NaN  NaN  NaN
     NaN      4.0  NaN    NaN    NaN       NaN  NaN  NaN    NaN  NaN  NaN  NaN
     NaN    NaN    NaN    NaN    NaN    …  NaN  NaN  NaN    NaN  NaN  NaN  NaN
       4.0  NaN    NaN    NaN    NaN       NaN  NaN  NaN    NaN  NaN  NaN  NaN
     NaN    NaN    NaN    NaN    NaN       NaN  NaN  NaN    NaN  NaN  NaN  NaN
     NaN    NaN    NaN    NaN    NaN       NaN  NaN  NaN    NaN  NaN  NaN  NaN
       4.0  NaN      1.0    2.0    1.0     NaN  NaN  NaN    NaN  NaN  NaN  NaN
     NaN    NaN    NaN      2.0  NaN    …  NaN  NaN  NaN    NaN  NaN  NaN  NaN
     NaN    NaN    NaN    NaN    NaN       NaN  NaN  NaN    NaN  NaN  NaN  NaN
     NaN    NaN    NaN    NaN    NaN       NaN  NaN  NaN    NaN  NaN  NaN  NaN
     NaN    NaN    NaN    NaN    NaN       NaN  NaN  NaN    NaN  NaN  NaN  NaN
       3.0  NaN    NaN    NaN    NaN       NaN  NaN  NaN    NaN  NaN  NaN  NaN
``` 


