make: x86_64-conda-linux-gnu-cc：命令未找到
conda install gxx_linux-64

bundle show minima
/home/jinyunfan/anaconda3/envs/jekyll-env/lib/ruby/gems/2.6.0/gems/minima-2.5.1

```{bash}
cp -r /home/jinyunfan/anaconda3/envs/jekyll-env/lib/ruby/gems/2.6.0/gems/minima-2.5.1/_layouts .
cp -r /home/jinyunfan/anaconda3/envs/jekyll-env/lib/ruby/gems/2.6.0/gems/minima-2.5.1/_includes .
cp -r /home/jinyunfan/anaconda3/envs/jekyll-env/lib/ruby/gems/2.6.0/gems/minima-2.5.1/assets .
```


- Reference
  - http://webdocs.cs.ualberta.ca/~zichen2/blog/coding/setup/2019/02/17/how-to-add-mathjax-support-to-jekyll.html
  - https://alan97.github.io/random/mathjax/

https://stackoverflow.com/questions/34422103/how-do-i-change-the-background-color-for-code-block-in-jekyll


- Install ruby 2.4
<https://blog.csdn.net/liguangxianbin/article/details/79454986>

`gem install jekyll bundler` gives
jekyll requires RubyGems version >= 2.7.0. Try 'gem update --system' to update RubyGems itself
- https://stackoverflow.com/questions/13626143/how-to-upgrade-rubygems
bundle update --bundler
https://github.com/jekyll/jekyll-seo-tag/blob/master/docs/installation.md


- 今天为啥总是没法push

```bash
# Run jekyll instance locally
bundle exec jekyll serve
```


```latex
\widetilde{} # A long ~
```

