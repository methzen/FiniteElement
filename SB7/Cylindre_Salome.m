function [Node,connect3D,g1,g2,g3]=Cylindre_Salome


Node=[
1 4.90600000000000e+00 0.00000000000000e+00 0.00000000000000e+00
2 5.00000000000000e+00 0.00000000000000e+00 0.00000000000000e+00
3 5.00000000000000e+00 0.00000000000000e+00 5.17500000000000e+00
4 4.90600000000000e+00 0.00000000000000e+00 5.17500000000000e+00
5 5.00000000000000e+00 0.00000000000000e+00 6.46875000000000e-01
6 5.00000000000000e+00 0.00000000000000e+00 1.29375000000000e+00
7 5.00000000000000e+00 0.00000000000000e+00 1.94062500000000e+00
8 5.00000000000000e+00 0.00000000000000e+00 2.58750000000000e+00
9 5.00000000000000e+00 0.00000000000000e+00 3.23437500000000e+00
10 5.00000000000000e+00 0.00000000000000e+00 3.88125000000000e+00
11 5.00000000000000e+00 0.00000000000000e+00 4.52812500000000e+00
12 4.90600000000000e+00 0.00000000000000e+00 4.52812500000000e+00
13 4.90600000000000e+00 0.00000000000000e+00 3.88125000000000e+00
14 4.90600000000000e+00 0.00000000000000e+00 3.23437500000000e+00
15 4.90600000000000e+00 0.00000000000000e+00 2.58750000000000e+00
16 4.90600000000000e+00 0.00000000000000e+00 1.94062500000000e+00
17 4.90600000000000e+00 0.00000000000000e+00 1.29375000000000e+00
18 4.90600000000000e+00 0.00000000000000e+00 6.46875000000000e-01
19 4.86402848989991e+00 6.40361499031573e-01 0.00000000000000e+00
20 4.73883210377417e+00 1.26976623527297e+00 0.00000000000000e+00
21 4.53255298650037e+00 1.87744491918313e+00 0.00000000000000e+00
22 4.24872063096646e+00 2.45300000000000e+00 0.00000000000000e+00
23 3.89219148746880e+00 2.98658357071678e+00 0.00000000000000e+00
24 3.46906586850120e+00 3.46906586850120e+00 0.00000000000000e+00
25 2.98658357071678e+00 3.89219148746880e+00 0.00000000000000e+00
26 2.45300000000000e+00 4.24872063096646e+00 0.00000000000000e+00
27 1.87744491918313e+00 4.53255298650037e+00 0.00000000000000e+00
28 1.26976623527297e+00 4.73883210377417e+00 0.00000000000000e+00
29 6.40361499031573e-01 4.86402848989991e+00 0.00000000000000e+00
30 2.22044604925031e-16 4.90600000000000e+00 0.00000000000000e+00
31 4.95722430686905e+00 6.52630961100258e-01 0.00000000000000e+00
32 4.82962913144534e+00 1.29409522551260e+00 0.00000000000000e+00
33 4.61939766255643e+00 1.91341716182545e+00 0.00000000000000e+00
34 4.33012701892219e+00 2.50000000000000e+00 0.00000000000000e+00
35 3.96676670145618e+00 3.04380714504360e+00 0.00000000000000e+00
36 3.53553390593274e+00 3.53553390593274e+00 0.00000000000000e+00
37 3.04380714504360e+00 3.96676670145617e+00 0.00000000000000e+00
38 2.50000000000000e+00 4.33012701892219e+00 0.00000000000000e+00
39 1.91341716182545e+00 4.61939766255643e+00 0.00000000000000e+00
40 1.29409522551260e+00 4.82962913144534e+00 0.00000000000000e+00
41 6.52630961100259e-01 4.95722430686905e+00 0.00000000000000e+00
42 7.77156117237610e-16 5.00000000000000e+00 0.00000000000000e+00
43 4.95722430686905e+00 6.52630961100258e-01 6.46875000000000e-01
44 4.82962913144534e+00 1.29409522551260e+00 6.46875000000000e-01
45 4.61939766255643e+00 1.91341716182545e+00 6.46875000000000e-01
46 4.33012701892219e+00 2.50000000000000e+00 6.46875000000000e-01
47 3.96676670145618e+00 3.04380714504360e+00 6.46875000000000e-01
48 3.53553390593274e+00 3.53553390593274e+00 6.46875000000000e-01
49 3.04380714504360e+00 3.96676670145617e+00 6.46875000000000e-01
50 2.50000000000000e+00 4.33012701892219e+00 6.46875000000000e-01
51 1.91341716182545e+00 4.61939766255643e+00 6.46875000000000e-01
52 1.29409522551260e+00 4.82962913144534e+00 6.46875000000000e-01
53 6.52630961100259e-01 4.95722430686905e+00 6.46875000000000e-01
54 7.77156117237610e-16 5.00000000000000e+00 6.46875000000000e-01
55 4.95722430686905e+00 6.52630961100258e-01 1.29375000000000e+00
56 4.82962913144534e+00 1.29409522551260e+00 1.29375000000000e+00
57 4.61939766255643e+00 1.91341716182545e+00 1.29375000000000e+00
58 4.33012701892219e+00 2.50000000000000e+00 1.29375000000000e+00
59 3.96676670145618e+00 3.04380714504360e+00 1.29375000000000e+00
60 3.53553390593274e+00 3.53553390593274e+00 1.29375000000000e+00
61 3.04380714504360e+00 3.96676670145617e+00 1.29375000000000e+00
62 2.50000000000000e+00 4.33012701892219e+00 1.29375000000000e+00
63 1.91341716182545e+00 4.61939766255643e+00 1.29375000000000e+00
64 1.29409522551260e+00 4.82962913144534e+00 1.29375000000000e+00
65 6.52630961100259e-01 4.95722430686905e+00 1.29375000000000e+00
66 7.77156117237610e-16 5.00000000000000e+00 1.29375000000000e+00
67 4.95722430686905e+00 6.52630961100258e-01 1.94062500000000e+00
68 4.82962913144534e+00 1.29409522551260e+00 1.94062500000000e+00
69 4.61939766255643e+00 1.91341716182545e+00 1.94062500000000e+00
70 4.33012701892219e+00 2.50000000000000e+00 1.94062500000000e+00
71 3.96676670145618e+00 3.04380714504360e+00 1.94062500000000e+00
72 3.53553390593274e+00 3.53553390593274e+00 1.94062500000000e+00
73 3.04380714504360e+00 3.96676670145617e+00 1.94062500000000e+00
74 2.50000000000000e+00 4.33012701892219e+00 1.94062500000000e+00
75 1.91341716182545e+00 4.61939766255643e+00 1.94062500000000e+00
76 1.29409522551260e+00 4.82962913144534e+00 1.94062500000000e+00
77 6.52630961100259e-01 4.95722430686905e+00 1.94062500000000e+00
78 7.77156117237610e-16 5.00000000000000e+00 1.94062500000000e+00
79 4.95722430686905e+00 6.52630961100258e-01 2.58750000000000e+00
80 4.82962913144534e+00 1.29409522551260e+00 2.58750000000000e+00
81 4.61939766255643e+00 1.91341716182545e+00 2.58750000000000e+00
82 4.33012701892219e+00 2.50000000000000e+00 2.58750000000000e+00
83 3.96676670145618e+00 3.04380714504360e+00 2.58750000000000e+00
84 3.53553390593274e+00 3.53553390593274e+00 2.58750000000000e+00
85 3.04380714504360e+00 3.96676670145617e+00 2.58750000000000e+00
86 2.50000000000000e+00 4.33012701892219e+00 2.58750000000000e+00
87 1.91341716182545e+00 4.61939766255643e+00 2.58750000000000e+00
88 1.29409522551260e+00 4.82962913144534e+00 2.58750000000000e+00
89 6.52630961100259e-01 4.95722430686905e+00 2.58750000000000e+00
90 7.77156117237610e-16 5.00000000000000e+00 2.58750000000000e+00
91 4.95722430686905e+00 6.52630961100258e-01 3.23437500000000e+00
92 4.82962913144534e+00 1.29409522551260e+00 3.23437500000000e+00
93 4.61939766255643e+00 1.91341716182545e+00 3.23437500000000e+00
94 4.33012701892219e+00 2.50000000000000e+00 3.23437500000000e+00
95 3.96676670145618e+00 3.04380714504360e+00 3.23437500000000e+00
96 3.53553390593274e+00 3.53553390593274e+00 3.23437500000000e+00
97 3.04380714504360e+00 3.96676670145617e+00 3.23437500000000e+00
98 2.50000000000000e+00 4.33012701892219e+00 3.23437500000000e+00
99 1.91341716182545e+00 4.61939766255643e+00 3.23437500000000e+00
100 1.29409522551260e+00 4.82962913144534e+00 3.23437500000000e+00
101 6.52630961100259e-01 4.95722430686905e+00 3.23437500000000e+00
102 7.77156117237610e-16 5.00000000000000e+00 3.23437500000000e+00
103 4.95722430686905e+00 6.52630961100258e-01 3.88125000000000e+00
104 4.82962913144534e+00 1.29409522551260e+00 3.88125000000000e+00
105 4.61939766255643e+00 1.91341716182545e+00 3.88125000000000e+00
106 4.33012701892219e+00 2.50000000000000e+00 3.88125000000000e+00
107 3.96676670145618e+00 3.04380714504360e+00 3.88125000000000e+00
108 3.53553390593274e+00 3.53553390593274e+00 3.88125000000000e+00
109 3.04380714504360e+00 3.96676670145617e+00 3.88125000000000e+00
110 2.50000000000000e+00 4.33012701892219e+00 3.88125000000000e+00
111 1.91341716182545e+00 4.61939766255643e+00 3.88125000000000e+00
112 1.29409522551260e+00 4.82962913144534e+00 3.88125000000000e+00
113 6.52630961100259e-01 4.95722430686905e+00 3.88125000000000e+00
114 7.77156117237610e-16 5.00000000000000e+00 3.88125000000000e+00
115 4.95722430686905e+00 6.52630961100258e-01 4.52812500000000e+00
116 4.82962913144534e+00 1.29409522551260e+00 4.52812500000000e+00
117 4.61939766255643e+00 1.91341716182545e+00 4.52812500000000e+00
118 4.33012701892219e+00 2.50000000000000e+00 4.52812500000000e+00
119 3.96676670145618e+00 3.04380714504360e+00 4.52812500000000e+00
120 3.53553390593274e+00 3.53553390593274e+00 4.52812500000000e+00
121 3.04380714504360e+00 3.96676670145617e+00 4.52812500000000e+00
122 2.50000000000000e+00 4.33012701892219e+00 4.52812500000000e+00
123 1.91341716182545e+00 4.61939766255643e+00 4.52812500000000e+00
124 1.29409522551260e+00 4.82962913144534e+00 4.52812500000000e+00
125 6.52630961100259e-01 4.95722430686905e+00 4.52812500000000e+00
126 7.77156117237610e-16 5.00000000000000e+00 4.52812500000000e+00
127 4.95722430686905e+00 6.52630961100258e-01 5.17500000000000e+00
128 4.82962913144534e+00 1.29409522551260e+00 5.17500000000000e+00
129 4.61939766255643e+00 1.91341716182545e+00 5.17500000000000e+00
130 4.33012701892219e+00 2.50000000000000e+00 5.17500000000000e+00
131 3.96676670145618e+00 3.04380714504360e+00 5.17500000000000e+00
132 3.53553390593274e+00 3.53553390593274e+00 5.17500000000000e+00
133 3.04380714504360e+00 3.96676670145617e+00 5.17500000000000e+00
134 2.50000000000000e+00 4.33012701892219e+00 5.17500000000000e+00
135 1.91341716182545e+00 4.61939766255643e+00 5.17500000000000e+00
136 1.29409522551260e+00 4.82962913144534e+00 5.17500000000000e+00
137 6.52630961100259e-01 4.95722430686905e+00 5.17500000000000e+00
138 7.77156117237610e-16 5.00000000000000e+00 5.17500000000000e+00
139 4.86402848989991e+00 6.40361499031573e-01 5.17500000000000e+00
140 4.73883210377417e+00 1.26976623527297e+00 5.17500000000000e+00
141 4.53255298650037e+00 1.87744491918313e+00 5.17500000000000e+00
142 4.24872063096646e+00 2.45300000000000e+00 5.17500000000000e+00
143 3.89219148746880e+00 2.98658357071678e+00 5.17500000000000e+00
144 3.46906586850120e+00 3.46906586850120e+00 5.17500000000000e+00
145 2.98658357071678e+00 3.89219148746880e+00 5.17500000000000e+00
146 2.45300000000000e+00 4.24872063096646e+00 5.17500000000000e+00
147 1.87744491918313e+00 4.53255298650037e+00 5.17500000000000e+00
148 1.26976623527297e+00 4.73883210377417e+00 5.17500000000000e+00
149 6.40361499031573e-01 4.86402848989991e+00 5.17500000000000e+00
150 2.22044604925031e-16 4.90600000000000e+00 5.17500000000000e+00
151 4.86402848989991e+00 6.40361499031573e-01 4.52812500000000e+00
152 4.73883210377417e+00 1.26976623527297e+00 4.52812500000000e+00
153 4.53255298650037e+00 1.87744491918313e+00 4.52812500000000e+00
154 4.24872063096646e+00 2.45300000000000e+00 4.52812500000000e+00
155 3.89219148746880e+00 2.98658357071678e+00 4.52812500000000e+00
156 3.46906586850120e+00 3.46906586850120e+00 4.52812500000000e+00
157 2.98658357071678e+00 3.89219148746880e+00 4.52812500000000e+00
158 2.45300000000000e+00 4.24872063096646e+00 4.52812500000000e+00
159 1.87744491918313e+00 4.53255298650037e+00 4.52812500000000e+00
160 1.26976623527297e+00 4.73883210377417e+00 4.52812500000000e+00
161 6.40361499031573e-01 4.86402848989991e+00 4.52812500000000e+00
162 2.22044604925031e-16 4.90600000000000e+00 4.52812500000000e+00
163 4.86402848989991e+00 6.40361499031573e-01 3.88125000000000e+00
164 4.73883210377417e+00 1.26976623527297e+00 3.88125000000000e+00
165 4.53255298650037e+00 1.87744491918313e+00 3.88125000000000e+00
166 4.24872063096646e+00 2.45300000000000e+00 3.88125000000000e+00
167 3.89219148746880e+00 2.98658357071678e+00 3.88125000000000e+00
168 3.46906586850120e+00 3.46906586850120e+00 3.88125000000000e+00
169 2.98658357071678e+00 3.89219148746880e+00 3.88125000000000e+00
170 2.45300000000000e+00 4.24872063096646e+00 3.88125000000000e+00
171 1.87744491918313e+00 4.53255298650037e+00 3.88125000000000e+00
172 1.26976623527297e+00 4.73883210377417e+00 3.88125000000000e+00
173 6.40361499031573e-01 4.86402848989991e+00 3.88125000000000e+00
174 2.22044604925031e-16 4.90600000000000e+00 3.88125000000000e+00
175 4.86402848989991e+00 6.40361499031573e-01 3.23437500000000e+00
176 4.73883210377417e+00 1.26976623527297e+00 3.23437500000000e+00
177 4.53255298650037e+00 1.87744491918313e+00 3.23437500000000e+00
178 4.24872063096646e+00 2.45300000000000e+00 3.23437500000000e+00
179 3.89219148746880e+00 2.98658357071678e+00 3.23437500000000e+00
180 3.46906586850120e+00 3.46906586850120e+00 3.23437500000000e+00
181 2.98658357071678e+00 3.89219148746880e+00 3.23437500000000e+00
182 2.45300000000000e+00 4.24872063096646e+00 3.23437500000000e+00
183 1.87744491918313e+00 4.53255298650037e+00 3.23437500000000e+00
184 1.26976623527297e+00 4.73883210377417e+00 3.23437500000000e+00
185 6.40361499031573e-01 4.86402848989991e+00 3.23437500000000e+00
186 2.22044604925031e-16 4.90600000000000e+00 3.23437500000000e+00
187 4.86402848989991e+00 6.40361499031573e-01 2.58750000000000e+00
188 4.73883210377417e+00 1.26976623527297e+00 2.58750000000000e+00
189 4.53255298650037e+00 1.87744491918313e+00 2.58750000000000e+00
190 4.24872063096646e+00 2.45300000000000e+00 2.58750000000000e+00
191 3.89219148746880e+00 2.98658357071678e+00 2.58750000000000e+00
192 3.46906586850120e+00 3.46906586850120e+00 2.58750000000000e+00
193 2.98658357071678e+00 3.89219148746880e+00 2.58750000000000e+00
194 2.45300000000000e+00 4.24872063096646e+00 2.58750000000000e+00
195 1.87744491918313e+00 4.53255298650037e+00 2.58750000000000e+00
196 1.26976623527297e+00 4.73883210377417e+00 2.58750000000000e+00
197 6.40361499031573e-01 4.86402848989991e+00 2.58750000000000e+00
198 2.22044604925031e-16 4.90600000000000e+00 2.58750000000000e+00
199 4.86402848989991e+00 6.40361499031573e-01 1.94062500000000e+00
200 4.73883210377417e+00 1.26976623527297e+00 1.94062500000000e+00
201 4.53255298650037e+00 1.87744491918313e+00 1.94062500000000e+00
202 4.24872063096646e+00 2.45300000000000e+00 1.94062500000000e+00
203 3.89219148746880e+00 2.98658357071678e+00 1.94062500000000e+00
204 3.46906586850120e+00 3.46906586850120e+00 1.94062500000000e+00
205 2.98658357071678e+00 3.89219148746880e+00 1.94062500000000e+00
206 2.45300000000000e+00 4.24872063096646e+00 1.94062500000000e+00
207 1.87744491918313e+00 4.53255298650037e+00 1.94062500000000e+00
208 1.26976623527297e+00 4.73883210377417e+00 1.94062500000000e+00
209 6.40361499031573e-01 4.86402848989991e+00 1.94062500000000e+00
210 2.22044604925031e-16 4.90600000000000e+00 1.94062500000000e+00
211 4.86402848989991e+00 6.40361499031573e-01 1.29375000000000e+00
212 4.73883210377417e+00 1.26976623527297e+00 1.29375000000000e+00
213 4.53255298650037e+00 1.87744491918313e+00 1.29375000000000e+00
214 4.24872063096646e+00 2.45300000000000e+00 1.29375000000000e+00
215 3.89219148746880e+00 2.98658357071678e+00 1.29375000000000e+00
216 3.46906586850120e+00 3.46906586850120e+00 1.29375000000000e+00
217 2.98658357071678e+00 3.89219148746880e+00 1.29375000000000e+00
218 2.45300000000000e+00 4.24872063096646e+00 1.29375000000000e+00
219 1.87744491918313e+00 4.53255298650037e+00 1.29375000000000e+00
220 1.26976623527297e+00 4.73883210377417e+00 1.29375000000000e+00
221 6.40361499031573e-01 4.86402848989991e+00 1.29375000000000e+00
222 2.22044604925031e-16 4.90600000000000e+00 1.29375000000000e+00
223 4.86402848989991e+00 6.40361499031573e-01 6.46875000000000e-01
224 4.73883210377417e+00 1.26976623527297e+00 6.46875000000000e-01
225 4.53255298650037e+00 1.87744491918313e+00 6.46875000000000e-01
226 4.24872063096646e+00 2.45300000000000e+00 6.46875000000000e-01
227 3.89219148746880e+00 2.98658357071678e+00 6.46875000000000e-01
228 3.46906586850120e+00 3.46906586850120e+00 6.46875000000000e-01
229 2.98658357071678e+00 3.89219148746880e+00 6.46875000000000e-01
230 2.45300000000000e+00 4.24872063096646e+00 6.46875000000000e-01
231 1.87744491918313e+00 4.53255298650037e+00 6.46875000000000e-01
232 1.26976623527297e+00 4.73883210377417e+00 6.46875000000000e-01
233 6.40361499031573e-01 4.86402848989991e+00 6.46875000000000e-01
234 2.22044604925031e-16 4.90600000000000e+00 6.46875000000000e-01];

Node=Node(:,2:4);

Con=[
243 308 1 2 5 18 19 31 43 223 
244 308 19 31 43 223 20 32 44 224 
245 308 20 32 44 224 21 33 45 225 
246 308 21 33 45 225 22 34 46 226 
247 308 22 34 46 226 23 35 47 227 
248 308 23 35 47 227 24 36 48 228 
249 308 24 36 48 228 25 37 49 229 
250 308 25 37 49 229 26 38 50 230 
251 308 26 38 50 230 27 39 51 231 
252 308 27 39 51 231 28 40 52 232 
253 308 28 40 52 232 29 41 53 233 
254 308 29 41 53 233 30 42 54 234 
255 308 18 5 6 17 223 43 55 211 
256 308 223 43 55 211 224 44 56 212 
257 308 224 44 56 212 225 45 57 213 
258 308 225 45 57 213 226 46 58 214 
259 308 226 46 58 214 227 47 59 215 
260 308 227 47 59 215 228 48 60 216 
261 308 228 48 60 216 229 49 61 217 
262 308 229 49 61 217 230 50 62 218 
263 308 230 50 62 218 231 51 63 219 
264 308 231 51 63 219 232 52 64 220 
265 308 232 52 64 220 233 53 65 221 
266 308 233 53 65 221 234 54 66 222 
267 308 17 6 7 16 211 55 67 199 
268 308 211 55 67 199 212 56 68 200 
269 308 212 56 68 200 213 57 69 201 
270 308 213 57 69 201 214 58 70 202 
271 308 214 58 70 202 215 59 71 203 
272 308 215 59 71 203 216 60 72 204 
273 308 216 60 72 204 217 61 73 205 
274 308 217 61 73 205 218 62 74 206 
275 308 218 62 74 206 219 63 75 207 
276 308 219 63 75 207 220 64 76 208 
277 308 220 64 76 208 221 65 77 209 
278 308 221 65 77 209 222 66 78 210 
279 308 16 7 8 15 199 67 79 187 
280 308 199 67 79 187 200 68 80 188 
281 308 200 68 80 188 201 69 81 189 
282 308 201 69 81 189 202 70 82 190 
283 308 202 70 82 190 203 71 83 191 
284 308 203 71 83 191 204 72 84 192 
285 308 204 72 84 192 205 73 85 193 
286 308 205 73 85 193 206 74 86 194 
287 308 206 74 86 194 207 75 87 195 
288 308 207 75 87 195 208 76 88 196 
289 308 208 76 88 196 209 77 89 197 
290 308 209 77 89 197 210 78 90 198 
291 308 15 8 9 14 187 79 91 175 
292 308 187 79 91 175 188 80 92 176 
293 308 188 80 92 176 189 81 93 177 
294 308 189 81 93 177 190 82 94 178 
295 308 190 82 94 178 191 83 95 179 
296 308 191 83 95 179 192 84 96 180 
297 308 192 84 96 180 193 85 97 181 
298 308 193 85 97 181 194 86 98 182 
299 308 194 86 98 182 195 87 99 183 
300 308 195 87 99 183 196 88 100 184 
301 308 196 88 100 184 197 89 101 185 
302 308 197 89 101 185 198 90 102 186 
303 308 14 9 10 13 175 91 103 163 
304 308 175 91 103 163 176 92 104 164 
305 308 176 92 104 164 177 93 105 165 
306 308 177 93 105 165 178 94 106 166 
307 308 178 94 106 166 179 95 107 167 
308 308 179 95 107 167 180 96 108 168 
309 308 180 96 108 168 181 97 109 169 
310 308 181 97 109 169 182 98 110 170 
311 308 182 98 110 170 183 99 111 171 
312 308 183 99 111 171 184 100 112 172 
313 308 184 100 112 172 185 101 113 173 
314 308 185 101 113 173 186 102 114 174 
315 308 13 10 11 12 163 103 115 151 
316 308 163 103 115 151 164 104 116 152 
317 308 164 104 116 152 165 105 117 153 
318 308 165 105 117 153 166 106 118 154 
319 308 166 106 118 154 167 107 119 155 
320 308 167 107 119 155 168 108 120 156 
321 308 168 108 120 156 169 109 121 157 
322 308 169 109 121 157 170 110 122 158 
323 308 170 110 122 158 171 111 123 159 
324 308 171 111 123 159 172 112 124 160 
325 308 172 112 124 160 173 113 125 161 
326 308 173 113 125 161 174 114 126 162 
327 308 12 11 3 4 151 115 127 139 
328 308 151 115 127 139 152 116 128 140 
329 308 152 116 128 140 153 117 129 141 
330 308 153 117 129 141 154 118 130 142 
331 308 154 118 130 142 155 119 131 143 
332 308 155 119 131 143 156 120 132 144 
333 308 156 120 132 144 157 121 133 145 
334 308 157 121 133 145 158 122 134 146 
335 308 158 122 134 146 159 123 135 147 
336 308 159 123 135 147 160 124 136 148 
337 308 160 124 136 148 161 125 137 149 
338 308 161 125 137 149 162 126 138 150 ];

%% corrections

Con=Con(:,3:10);

iter=size(Con,1);
f=zeros(1,8);

for i=1:iter
    g=Con(i,:);
    f(1)=g(1); f(2)=g(4); f(3)=g(8); f(4)=g(5);
    f(5)=g(2); f(6)=g(3); f(7)=g(7); f(8)=g(6);
    Con(i,:)=f;   
    
end

a=Con(:,1:4);
Con(:,1:4)=Con(:,5:8);
Con(:,5:8)=a;
%%

connect3D=Con;

g1=[30 % x
    42
    54
    66
    78
    90
   102
   114
   126
   138
   150
   162
   174
   186
   198
   210
   222
   234];
g2=[ 1 % y
     2
     3
     4
     5
     6
     7
     8
     9
    10
    11
    12
    13
    14
    15
    16
    17
    18];

g3=[ 1 %z
     2
    19
    20
    21
    22
    23
    24
    25
    26
    27
    28
    29
    30
    31
    32
    33
    34
    35
    36
    37
    38
    39
    40
    41
    42];
end