---
layout: mathPage
title: Hirschberg Algorithm
---

The implementation main file: [Hirschberg.py](https://github.com/M-Abdallah/algorithms-sbe631/blob/master/Hirschberg.py){:target="_blank"}

This is similar to the Needleman--Wuesch score grid calculation based on  match and and mismatch; the  only difference is that it does not do it for the whole sequence.

``` python
def scoreCalc( seq1, seq2 ):
    # This calculate the of the scoreGrid for those sequences seq1 and seq2
    scoreGrid = np.zeros((seq1.size+1, seq2.size+1), dtype=np.int16)
    # Fill the score grid(matching grid)
    for i in range(seq1.size+1):
        scoreGrid[i, 0] = i*gapPen
    for j in range(seq2.size+1):
        scoreGrid[0, j] = j*gapPen
    # Filling using the score matrix
    for i in range(1, seq1.size+1):
        for j in range(1, seq2.size+1):
            scoreGrid[i, j] = max(scoreGrid[i-1, j-1] + getScore(seq1, seq2, i, j),
                                  scoreGrid[i-1, j] + gapPen,
                                  scoreGrid[i, j-1] + gapPen)
    
    return scoreGrid
```

This just return the last row in the array to be used later in the divide part.
```python
def NWscore(seq1, seq2):
    # This calculate the last row of the of the scoreGrid for those sequences seq1 and seq2
    scoreGrid = scoreCalc( seq1, seq2 )
    return scoreGrid[seq1.size, :]
```

The following is where all the fun happens 

- This part adds the gaps in either the first sequence and the second sequence
```python
def traceback(seqX, seqY):
    # Translate the score(matching) grid to matched pairs or gaps
    alignedPair = []
    i = seqX.size
    j = seqY.size
    # scoreGrid = scoreCalc( seqX, seqY )
    if j == 0:
        # This is search for direct mismatch and add gaps for the second sequence
        while i > 0:
            alignedPair.append(getAlignedPair(seqX, seqY, i, 0))
            i -= 1
    elif i == 0:
        # This is search for direct mismatch and add gaps for the first sequence
        while j > 0:
            alignedPair.append( getAlignedPair(seqX, seqY,0, j))
            j -= 1
```

- This is just normal Needleman-Wunsch Algorithm, but on a smaller array.
```python
    elif i == 1 or j == 1:
        # This is just normal Needleman-Wunsch Algorithm 
        scoreGrid = scoreCalc( seqX, seqY )
        while i > 0 or j > 0:
            if i > 0 and j > 0 and scoreGrid[i-1, j-1] + getScore( seqX,seqY, i, j) == scoreGrid[i, j]:
                alignedPair.append(getAlignedPair( seqX,seqY, i, j))
                i -= 1
                j -= 1
            elif i > 0 and scoreGrid[i-1, j] + gapPen == scoreGrid[i, j]:
                alignedPair.append( getAlignedPair( seqX,seqY, i, 0))
                i -= 1
            else:
                alignedPair.append(getAlignedPair(seqX, seqY, 0, j))
                j -= 1
```

- We divide the data in this part.
    1. The first sequence is divided normally, by dividing the array by 2.
    2. The second sequence is divided based the location of the maximum value of the last row of the score grid for these sequences.
After that we call the trace back function to align the two parts.

``` python
    else:
        # Here is were all the magic happens, and by magic I mean Divide&Conquer happens
        # Where the First sequence is divided
        iMid = int(i/2)
        scoreL = NWscore( seqX[0:iMid], seqY )
        scoreR = NWscore( seqX[iMid:i][::-1], seqY[::-1] )
        sum = np.add(scoreL, scoreR[::-1] )
        
        jMid = np.argmax( sum )
        alignedPair+= traceback ( seqX[0:iMid], seqY[ 0: jMid ] )
        alignedPair+= traceback ( seqX[iMid:i], seqY[ jMid: j ] )  
    
    print(alignedPair)
    return alignedPair
```