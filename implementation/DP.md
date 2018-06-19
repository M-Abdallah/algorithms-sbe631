---
layout: mathPage
title: Needleman-Wuesch Algorithm
---

The implementation main file: [Needleman-Wunsch_Alg.py](https://github.com/M-Abdallah/algorithms-sbe631/blob/master/Needleman-Wunsch_Alg.py){:target="_blank"}

The scoring matrix AKA [Similarity Matrix](https://en.wikipedia.org/wiki/Similarity_measure). Basically, it score the mismatches as a -ve values and the matches ans +ve values.


``` python
gapPen = -1

# scoring = np.array([[10,-3,-1,-4],
#                     [-3,9,-5,-0],
#                     [-1,-5,7,-3],
#                     [-4,-0,-3,8]])

scoring = np.array([[3,-3,-3,-3],
                    [-3,3,-3,-3],
                    [-3,-3,3,-3],
                    [-3,-3,-3,3]])
```

gapPen ← MismatchScore

$$
scoreGrid(0,j) = gapPen ∗ j
$$

$$
scoreGrid(i,0) = gapPen ∗ i
$$

Recursion, based on the principle of optimality:

$$
scoreGrid(i,j) = max( scoreGrid(i-1,j-1) + scoreMatrix(A_i,B_j), scoreGrid(i,j-1) + gapPen, scoreGrid(i-1,j) + gapPen )
$$

``` python
    def __computeArray__(self):
        # Fill the score grid(matching grid)
        for i in range(self.seq1.size+1):
            self.scoreGrid[i,0] = i*gapPen
        for j in range(self.seq2.size+1):
            self.scoreGrid[0,j] = j*gapPen
        # Filling using the score matrix 
        for i in range(1, self.seq1.size+1):
            for j in range(1, self.seq2.size+1):
                self.scoreGrid[i,j] = max(  self.scoreGrid[i-1, j-1] + self.__getScore__(i, j),
                                    self.scoreGrid[i-1, j] + gapPen,
                                    self.scoreGrid[i, j-1] + gapPen)
```

Using ```python
A, C, G, T = 0, 1, 2, 3```
we can find the score based on the elements of the two sequences

``` python
    def __getScore__(self, i, j):
        # we check it fin its score based on the matching(score) matrix we defined first and the i-1&j-1 
        return scoring[self.seq1[i-1], self.seq2[j-1]]
```

```python    
    def __getAlignedPair__(self, i, j):
        n1 = int_to_char[self.seq1[i-1]] if i>0 else '_'
        n2 = int_to_char[self.seq2[j-1]] if j>0 else '_'
        return (n1, n2)
```

``` python
    def __traceback__(self):
        # Translate the score(matching) grid to matched pairs or gaps
        alignedPair= []
        i = self.seq1.size
        j = self.seq2.size
        while i >0 or j>0:
            if i >0 and j>0 and self.scoreGrid[i-1, j-1] + self.__getScore__(i, j) == self.scoreGrid[i,j]:
                alignedPair.append(self.__getAlignedPair__(i, j))
                i -= 1
                j -= 1
            elif i >0 and self.scoreGrid[i-1, j] + gapPen == self.scoreGrid[i,j]:
                alignedPair.append(self.__getAlignedPair__(i, 0))
                i -= 1
            else:
                alignedPair.append(self.__getAlignedPair__(0, j))
                j -= 1
        alignedPair.reverse()
        return alignedPair  
```