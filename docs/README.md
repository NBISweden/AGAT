### Mkdocs

#### Welcome to MkDocs

For full documentation visit [mkdocs.org](https://www.mkdocs.org).  
For full documentation about the [material mkdocs theme](https://squidfunk.github.io/mkdocs-material/).

#### Installation

##### Manual

As prerequisite you need python >=3.8 and pip.  

Install Mkdocs:

`pip install mkdocs`

For the theme:  
`pip install mkdocs-material`

For the extensions:  
`pip install pymdown-extensions`

For the plugins:  
`pip install mkdocs-minify-plugin`  
`pip install mkdocs-macros-plugin`
`pip install mkdocs-embed-external-markdown`

##### Conda

Clone the repository and move in it.  
Then install all dependencies using conda and the `conda_env.yml` shipped with this repo:

```
conda env create -f conda_env.yml
```

Activate the environment and you are good:

```
conda activate education
```

#### Testing and building the website


* `mkdocs serve` - Start the live-reloading docs server, to test the site locally (http://127.0.0.1:8000/).
* `mkdocs gh-deploy` - Deploys the site on github pages.

* `mkdocs build` - Build the documentation site.
* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs -h` - Print help message and exit.
