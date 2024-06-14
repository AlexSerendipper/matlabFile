# Latex
\usepackage{amsmath,amsfonts}
\usepackage{algorithmic}
\usepackage{algorithm}
\usepackage{array}
\usepackage[caption=false,font=normalsize,labelfont=sf,textfont=sf]{subfig}
\usepackage{textcomp}
\usepackage{stfloats}
\usepackage{url}
\usepackage{verbatim}
\usepackage{graphicx}
\usepackage{cite}
\usepackage{textcomp}
\usepackage{upgreek}
\usepackage{multirow}
\usepackage{multicol}

## 0. 作者域

## 1. 公式域

#### 单公式单编号
``` latex
\begin{equation}
    \label{eq9}
       ...
       ...
\end{equation}
```

#### 多行公式多编号
- 使用\nonumber可以实现某行没有编号
- 注意 \begin{align} 与 \begin{equation}不能同时使用
- 使用\begin{align*}，可以使用多行公式无编号
- align，split，cases 的公式，都可以使用&进行对齐
<img src="image/latex1.png" style="zoom:60">
``` latex
\begin{align}
    \hat{C}_t &= tanh(W_C\odot[h_{t-1},x_t]+b_C) \\
    i_t&=\sigma(W_i\odot[h_{t-1},x_t]+b_i)  \nonumber \\ 
    f_t&=\sigma(W_f\odot[h_{t-1},x_t]+b_f) \\
    C_t&=f_t * C_{t-1}+i_t * \hat{C}_t \\
    o_t&=\sigma(W_o\odot[h_{t-1},x_t]+b_o) \\
    h_t&=o_t * tanh(C_t)
\end{align}
```

#### 多行公式单编号，编号居中
<img src="image/latex2.png" style="zoom:60">
``` latex
\begin{equation}
    \begin{aligned}
        &\hat{C}_t = tanh(W_C\odot[h_{t-1},x_t]+b_C) \\
        &i_t=\sigma(W_i\odot[h_{t-1},x_t]+b_i) \\ 
        &f_t=\sigma(W_f\odot[h_{t-1},x_t]+b_f) \\
        &C_t=f_t * C_{t-1}+i_t * \hat{C}_t \\
        &o_t=\sigma(W_o\odot[h_{t-1},x_t]+b_o) \\
        &h_t=o_t * tanh(C_t)
    \end{aligned}
\end{equation}
```
#### 一行公式分多行写 
<img src="image/latex3.png" style="zoom:60">
- align，split，cases 的公式，都可以使用&进行对齐
``` latex
\begin{equation}
	\begin{split}
	\cos 2x &= \cos^2 x - \sin^2 x\\
	        &= 2\cos^2 x - 1
	\end{split}
\end{equation}
```

#### 按情况（case）划分 
- 该方式和直接使用大括号类似
<img src="image/latex4.png" style="zoom:60">
``` latex
\begin{equation}
	D(x) = \begin{cases}
	      1, & if \ x \in Q \\
	      0, & if \ x \in R	
		   \end{cases}
\end{equation}
```

#### 大括号单编号
- \left\{ 表示添加左大括号，\right. 表示不加右括号结束
<img src="image/latex4.png" style="zoom:60">
``` latex
\begin{align}
    \left\{
        \begin{aligned}
        x&=eq1\\
        y&=eq2+1
        \end{aligned}
    \right.
\end{align}
````
#### 大括号多编号
<img src="image/latex5.png" style="zoom:60">
``` latex
\begin{numcases}{}
	x_1&=eq1 \label{eq1} \\
	x_2+1&=eq2 \label{eq2}
\end{numcases}
````
#### 调整公式内字体大小（替换small即可）
``` latex
\begin{small}
    \begin{equation}
        ...
    \end{equation}
\end{small}

七号 　　5.25pt 　　 1.845mm　　　　\tiny
六号 　　7.875pt 　　 2.768mm　　　　\scriptsize
小五号 　9pt 　　　　 3.163mm　　　　\footnotesize
五号 　　10.5pt 　　 3.69mm　　　　 \small
小四号 　12pt 　　　　4.2175mm　　　 \normalsize
四号 　　13.75pt 　　 4.83mm　　　　 \large
三号 　　15.75pt 　　 5.53mm　　　　 \Large
二号 　　21pt 　　　　7.38mm        \LARGE
一号 　　27.5pt 　　 9.48mm　　　　 \huge
小初号 　36pt 　　　　12.65mm　　　　\Huge
```
## 2. 图片域
#### 普通插入图片
- !h表示在当前处插入，!t表示在当前页的顶部插入图片，!b表示在当前页的底部插入图片
- centering意味着将图片居中
- 设置图片大小
   - scale=0.6，表示按照原图片长宽比的0.6倍进行缩放。。。
   - width=宽度,height=高度,keepaspectratio。。。可以指定长宽，并用keepaspectratio锁定纵横比
   - width=0.5\textwidth，图片的宽度被设置为了文本宽度的一半。你也可以用其他单位，比如cm或in。
   - width=\linewidth表示图片大小自适应文档宽度
- label命令添加了一个标签，这样你就可以在文中使用\ref{fig:example}来引用这个图片。
- \begin{figure\*}，\end{figure\*},通过在figure后加*号，可以使用图片横跨两栏，常用于在双栏模板下，图片放置两栏

``` latex
\begin{figure}[!h]
    \centering
    \includegraphics[scale=0.3]{Figs/fig5.png}
    \caption{图片附文}
    \label{fig_5}
\end{figure}
```
#### 并列插入图片（解决不了，有点问题一直报错）

## 3. 表格域
- Table：表格的所有元素
- Tabular：仅指表格的那部分，不包括表格的文字说明等
   - tabular处就可以指定表格排版，如{l|c}，l表示左对齐，c表示居中对齐，r表示右对齐，| 表示两列之间的分隔符，当不加|时， 这两列之间就不会有竖线。
- 表格中
   - `&`：表示一行中两个单元格之间的分隔符。比如，这个有两列，则每行有一个&符。如果有n列，则每行应该有 n-1 个 & 符
   - `\\`：表示换行符，表示该行结束，换下一行
   - `\hline`：horizontal line，表示在该行下面应该增加一条水平线。
   - `\cline{x-y}`:合并单元格时使用，在改行下的第x列到y列 增加一条水平线
   - `\multirow{行数}{*}{文本}`，行合并单元格
   - `\multicolumn{列数}{c}{文本}`，列合并单元格。此处可以指定合并后列的分隔符，如`c|`
   

<img src="image/latex15.png" style="zoom:60">
``` latex
\begin{table}[!h]
\caption{}
\label{tab:table1}
\centering
\begin{tabular}{|c|c|c|c|c|c|c|c|}
    \hline
    \multirow{2}{*}{$C$} & \multicolumn{7}{c|}{MS error based on different lengths of window function(nm)} \\ \cline{2-8}
    ~ & 64pt & 128pt & 256pt & 512pt & 1024pt & 2048pt & 4096pt \\ \hline
    0.1 & 14.93 & 11.04 & 5.47 & 3.41 & 6.02 & 532.25 & 165.40 \\ \hline
    0.5 & 15.68 & 12.33 & 7.00 & 5.91 & 7.90 & 459.98 & 207.43 \\ \hline
    0.9 & 17.35 & 14.64 & 9.39 & 7.57 & 12.37 & 737.04 & 481.09 \\ \hline
    1.2 & 18.86 & 16.69 & 12.14 & 8.87 & 14.98 & 581.63 & 2412.20 \\ \hline
    1.5 & 22.40 & 20.97 & 17.56 & 14.48 & 21.43 & 2276.81 & 8271.94 \\ \hline
\end{tabular}
\end{table}
```

## 4. 算法域

<img src="image/latex14.png" style="zoom:60">

``` latex
\begin{algorithm}[H]
\caption{Algorithm 1 $L(t)=TSPM(P_0(t),l,r)$ }
\label{alg:alg1}
...
...
\label{alg1}
\end{algorithm}
```

## 5.参考文献域
1. .bib文件中，使用逗号对作者进行分割是不合法的，只能通过`and`以及`\&`的方式对作者名进行分割，并且`and`是会被翻译成`,`的
2. 参考文献格式
- 在文中引入参考文献的格式是 \cite{bib1, bib2, bib3}，经过编译后得到的的形式是[1][2][3]，然而这并不是我们想要的格式，我们想要的格式是[1-3], 因此请参考第二种形式。
- 在Latex源文件 .tex文档的顶部引入包：\usepackage[numbers,sort&compress]{natbib}
在文中引入参考文献：\cite{bib1, bib2, bib3} ，经过编译后得到了我们想要的参考文献格式[1-3]
3. Latex参考文献中大写字母编译后自动变成了小写，可以将要保持原大小的内容置于{}内，这样就能保证不会被翻译成小写内容
4. author中有三种引用方式
（1）使用逗号时：姓氏在前，名在后，并用逗号隔开，and连接两个人名
- author = {Mohri, Mehryar and Rostamizadeh, Afshin and Talwalkar, Ameet},
（2）不使用逗号时：按照西文正常写法，名在前，姓在后，and连接两个人名
- author = {Mehryar Mohri and Afshin Rostamizadeh and Ameet Talwalkar}, 
（3）手工硬编码缩写：按照西文正常写法，名在前，姓在后，注意姓名中间有空格，and连接两个人名
- author = {M. Mohri and A. Rostamizadeh and A. Talwalkar},
5. 使用and others编写，会自动翻译成et al.
- author={Y. {Jiang} and X. {Gong} and D. {Liu} and others}, 
## 3. 公式中常用，$ $

| Option | Description |
| ------:| -----------:|
| α | \alpha |
| δ | \delta |
|△| \Delta |
| Λ | \Lambda |
| ω | \omega |
| σ | \sigma |
| χ | \chi |
| τ | \tau |
| ψ | \psi |
| Φ | \phi |
| υ | \upsilon |
| ν | \nu |
| u | \upmu（微米的第一个字母） |
| ε | \epsilon |
| ∈ | \in |
| 长→ | \longrightarrow |
| <img src="image/latex12.png" style="zoom:0"> | \hat{P_0} |
| · | \cdot(点乘) |
| x | \times(叉乘) |
| ≈ | \approx |
| ∼ | \sim（在公式中直接输入~号是代表空格的意思） |
| − | 注意，公式中或者文章中，都是加载不了这个减号的，暂不知怎么打出来这个。要改成- | | ^ | 上标 |
| _ | 下标（使用花括号进行连续的下标） |


## 3. 其他常用
| Option | Description |
| ------:| -----------:|
| 求和符号，下标为1，上标为无穷 | \sum_{n=1}^{\infty} |
| 分数 | \frac{分子}{分母} |
| <img src="image/latex16.png" style="zoom:0"> | ^{18}O |
| <img src="image/latex6.png" style="zoom:0"> | \left( \frac{a}{b} \right) |
| <img src="image/latex7.png" style="zoom:0"> | \left[ \frac{a}{b} \right] |
| <img src="image/latex8.png" style="zoom:0"> | \left\{ \frac{a}{b} \right\} |
| <img src="image/latex9.png" style="zoom:0"> | \left|
| <img src="image/latex10.png" style="zoom:0"> | \left \| \frac{a}{b} \right \| |
| <img src="image/latex11.png" style="zoom:0"> | \left \{ \frac{a}{b} \right . |
| 两个符号叠加<img src="image/latex13.png" style="zoom:0"> | f_n(x) \stackrel{*}{\approx} 1 |
| 字体加粗（主要用于公式内） | \usepackage{bm}。。bm{}=\mathbf{}|
| 字体加粗（文本域使用） | \textbf{}|
|  <img src="image/latex17.png" style="zoom:0"> | \sum_{i=0}^n|

