@def title = "Converting latex/pdf to html/markdown  for Jekyll-based blogs"
@def published ="January 6, 2023"
@def tags =["programming"]

# Converting latex/pdf to html/markdown for Jekyll-based blogs

**Shuvomoy Das Gupta**

*January 6, 2023*

In this blog, we will discuss how to convert a `.tex` or `.pdf` file to `.html` and then use that `.html` in a markdown blog post.

#### `.tex` to `.md` 

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

#### `.pdf` to `.html` using `pdf2htmlEX` 

#### On Windows 

*  download `pdf2htmlEX` zip file from the link [https://shuvomoy.github.io/blogs/assets/pdf2htmlEX/pdf2htmlEX-win32-0.14.6-with-poppler-data.zip](https://shuvomoy.github.io/blogs/assets/pdf2htmlEX/pdf2htmlEX-win32-0.14.6-with-poppler-data.zip)
* extract it to a location of your choice, eg, `C:\Program Files (x86)\ `
* now copy the pdf (say `filename.pdf`) file to the same folder `C:\Program Files (x86)\pdf2htmlEX-win32-0.14.6-with-poppler-data` which contains the `pdf2htmlEX.exe`
* in that folder, right-click to open a terminal (it will open `windows powershell`) and then run ` .\pdf2htmlEX --zoom 1.75 --process-outline=0 ".\filename.pdf"`
* (in case pdf2htmlEX cannot render the file, the ultimate hand is:  ` .\pdf2htmlEX --zoom 1.75 --fallback 1 --process-outline=0 ".\filename.pdf"`, which will ensure maximum possible accuracy and compatibility at the expense of a larger file size) 
* which will create a `filename.html` in the same folder
* copy the file in your `posts` folder, which will have url such as: `http://localhost:8000/posts/filename/`

#### On Linux

For converting pdf to html in a linux based OS (or [Linux on Windows with WSL](https://learn.microsoft.com/en-us/windows/wsl/install)), do the following steps:

* install `ttfautohint` by inputting the following in terminal `sudo apt install ttfautohint`

* install `pdf2htmlEX` from the link [https://shuvomoy.github.io/blogs/assets/pdf2htmlEX/pdf2htmlEX-0.18.8.rc1-master-20200630-Ubuntu-focal-x86_64.deb] (download the deb and then install it via package manager or by inputting `sudo dpkg -i pdf2htmlEX-0.18.8.rc1-master-20200630-Ubuntu-focal-x86_64.deb.deb`)
* go to the folder containing the pdf file by typing `cd DIR_NAME` in terminal
* convert the file into html format by typing 
 `pdf2htmlEX --zoom 1.75 --external-hint-tool=ttfautohint --process-outline=0 "filename.pdf"`, 
* which will create the `filename.html` file
* copy the file in your `posts` folder, which will have url such as: `http://localhost:8000/posts/filename/`
