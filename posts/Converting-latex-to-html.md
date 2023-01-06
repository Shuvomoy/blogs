@def title = "Converting latex to html  for Jekyll-based blogs"
@def published ="January 6, 2023"
@def tags =["programming"]

# Converting latex to html to markdown for Jekyll-based blogs

**Shuvomoy Das Gupta**

*January 6, 2023*

In this blog, we will discuss how to convert a `.tex` file to `.html` and then use that `.html` in a markdown blog post.

Suppose, name of the  `.tex` file is `test.tex`. First, change the directory to the one where this file resides. Type the following in the terminal:

```
cd "path_of_the_folder_containing_test.tex"
make4ht test.tex "mathml,mathjax"
```

Here, `make4ht ` is a build system that comes with [texlive](https://www.tug.org/texlive/), so please install the latter if you do not have it already. The last command will create a file named `test.html` in the same folder. Now open the `test.html` file in a text editor, and copy the entire content of the file. Now open the markdown file, which would be used to create the blog post, and then paste the copied content between the two `~~~`s as follows.

```
~~~
paste the content of test.html
~~~
```

That's it, the markdown file can be used as a blog post in Jekyll based website. 
