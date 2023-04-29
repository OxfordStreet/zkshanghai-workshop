# 第2课 课后作业

## circom编程练习：https://zkshanghai.xyz/notes/exercise2.html

### 转换为bit位 Num2Bits

* 参数：`nBits`
* 输入信号：`in`
* 输出信号：`b[nBits]`

输出信号应该是长度为 `nBits`的位数组，相当于 `in`的二进制表示。 `b[0]` 是最低有效位。

[解决方案](https://github.com/iden3/circomlib/blob/master/circuits/bitify.circom#L25)

```
pragma circom 2.0.0;

template Num2Bits(n) {
    signal input in;
    signal output out[n];

    var e2 = 1;
    for (var i = 0; i < n; i++) {
        out[i] <-- (in >> i) & 1;
        enforce out[i] == 0 || out[i] == 1;
    }

    enforce in === multiExp(out, [e2...e2*(n-1)]);
}

component main {
    signal private a;
    signal public out[8];

    // Connect input to template
    Num2Bits(8) num2bits = Num2Bits(8)();
    connect(a, num2bits.in);
    connect(num2bits.out, out);

    // Set default value for a
    a <== 100;
}

```

### 判零 IsZero

* 参数：无
* 输入信号：`in`
* 输出信号：`out`

要求：如果 `in`为零，`out`应为 `1`。 如果 `in`不为零，`out`应为 `0`。 这个有点棘手！

[解决方案](https://github.com/iden3/circomlib/blob/master/circuits/comparators.circom#L24)

```
pragma circom 2.1.4;

template IsZero() {
    signal input in;
    signal output out;

    signal inv;

    inv <-- in != 0 ? 1 / in : 0;

    out <== -in * inv + 1;
    in * out === 0;
}

component main {
    // No inputs necessary
    signal output out;
  
    IsZero isZero = IsZero();
  
    // Connect the output of IsZero to the main output
    isZero.out -> out;
} = main;

```

### 相等 IsEqual

* 参数：无
* 输入信号：`in[2]`
* 输出信号：`out`

要求：如果 `in[0]` 等于 `in[1]`，则 `out` 应为 `1`。 否则，`out` 应该是 `0`。

[解决方案](https://github.com/iden3/circomlib/blob/master/circuits/comparators.circom#L37)

```
pragma circom 2.1.4;

template IsEqual() {
    signal input in[2];
    signal output out;

    out <== in[0] == in[1] ? 1 : 0;
}

component main {
    // No inputs necessary
    signal input in[2];
    signal output out;
  
    IsEqual isEqual = IsEqual();
  
    // Connect the inputs of main to the inputs of IsEqual
    in[0] -> isEqual.in[0];
    in[1] -> isEqual.in[1];
  
    // Connect the output of IsEqual to the main output
    isEqual.out -> out;
} = main;

```

### 选择器 Selector

* 参数：`nChoices`
* 输入信号：`in[nChoices]`, `index`
* 输出：`out`

要求：输出 `out`应该等于 `in[index]`。 如果 `index` 越界（不在 [0, nChoices) 中），`out` 应该是 `0`。

[解决方案](https://github.com/darkforest-eth/circuits/blob/master/perlin/QuinSelector.circom)

```
pragma circom 2.1.4;

include "CalculateTotal.circom"

template QuinSelector(choices) {
    signal input in[choices];
    signal input index;
    signal output out;
  
    // Ensure that index < choices
    component lessThan = LessThan(4);
    lessThan.in[0] <== index;
    lessThan.in[1] <== choices;
    lessThan.out === 1;

    component calcTotal = CalculateTotal(choices);
    component eqs[choices];

    // For each item, check whether its index equals the input index.
    for (var i = 0; i < choices; i ++) {
        eqs[i] = IsEqual();
        eqs[i].in[0] <== i;
        eqs[i].in[1] <== index;

        // eqs[i].out is 1 if the index matches. As such, at most one input to
        // calcTotal is not 0.
        calcTotal.in[i] <== eqs[i].out * in[i];
    }

    // Returns 0 + 0 + 0 + item
    out <== calcTotal.out;
}

component main {} = QuinSelector();


/* INPUT = {
    "a": "5",
    "b": "77"
} */
```

错误信息有待debug：

```
STDERR: 
error[P1000]: UnrecognizedToken { token: (55, Token(71, "template"), 63), expected: ["\")\"", "\",\"", "\";\""] }
  ┌─ "main.circom":5:1
  │
5 │ template QuinSelector(choices) {
  │ ^^^^^^^^ Invalid syntax

previous errors were found
STDOUT: 
Compiled in 0.98s
LOG: 
hash 6008246173323011098915936938805752727781568490715388424063708882447636047656
OUTPUT: 
c = 385
ARTIFACTS: 
Finished in 1.10s
main.wasm (1027.06KB)
main.js (9.18KB)
main.wtns (7.88KB)
main.r1cs (110.08KB)
main.sym (37.71KB)
```

### 判负 IsNegative

注意：信号是模 p（Babyjubjub 素数）的残基，并且没有 `负`数模 p 的自然概念。 但是，很明显，当我们将 `p-1`视为 `-1`时，模运算类似于整数运算。 所以我们定义一个约定：`取负` 按照惯例认为是 (p/2, p-1] 中的余数，非负是 [0, p/2) 中的任意数

* 参数：无
* 输入信号：`in`
* 输出信号：`out`

要求：如果根据我们的约定，`in` 为负数，则 `out` 应为 `1`。 否则，`out` 应该是 `0`。 您可以自由使用[CompConstant circuit](https://github.com/iden3/circomlib/blob/master/circuits/compconstant.circom)，它有一个常量参数 `ct`，如果 `in`（二进制数组）在解释为整数时严格大于 `ct` 则输出 `1` ，否则为 `0`。

[解决方案](https://github.com/iden3/circomlib/blob/master/circuits/sign.circom#L23)

```
pragma circom 2.1.4;

// Define the p value for Babyjubjub
define P = 21888242871839275222246405745257275088548364400416034343698204186575808495617n;

template IsNegative() {
    signal input in[254];
    signal output out;
  
    CompConstant compare = CompConstant(ct=to_bits(P/2));
  
    // Negate the input signal by subtracting it from P
    signal[] negIn = mulmod(in, to_bits(P-1), to_bits(P));
  
    // Determine if the negative input is less than P/2
    signal lessThanHalf = compare(in=negIn);
  
    out <== lessThanHalf;
}

component main {
    signal input in[254];
    signal output out;
  
    IsNegative isNegative = IsNegative();
  
    in -> isNegative.in;
    isNegative.out -> out;
} = main;
```

* **理解检查** ：为什么我们不能只使用 LessThan 或上一个练习中的比较器电路之一？

  因为在 Babyjubjub 模中 `负`数的概念和整数不太一样，并且没有简单的方法处理它们，所以我们使用了 `CompConstant` 组件来判断是否小于实际上是大于相应的负数。这种方式可以避免直接比较负数本身的问题。此外，在这个例子中，我们需要额外的乘法模运算操作，以将输入信号取负。因此，使用 `CompConstant` 可以简化电路并更好地解决问题。

### 少于 LessThan

* 参数：无
* 输入信号：`in[2]`。 假设提前知道这些最多 2252−1。
* 输出信号：`out`

要求：如果 `in[0]` 严格小于 `in[1]`，则 `out` 应为 `1`。 否则，`out` 应该是 `0`。

* **扩展 1** ：如果您知道输入信号最多为 2k−1(k≤252)，您如何减少该电路所需的约束总数？ 编写一个在 `k`中参数化的电路版本。
* **扩展 2** ：编写 LessEqThan（测试 in[0] 是否 ≤ in[1]）、GreaterThan 和 GreaterEqThan

[解决方案（扩展1）](https://github.com/iden3/circomlib/blob/master/circuits/comparators.circom#L89)

```
include "utils.circom";

template LessThan() {
  signal in[2];
  signal out;

  // 将输入信号拆分成单独的比特。
  signal in_bit0 = bits(in[0], 32);
  signal in_bit1 = bits(in[1], 32);

  // 比较输入信号的每个比特。
  component cmps[32] = [0; 32];
  for (i = 0; i < 32; i++) {
    cmps[i] = GT();
    cmps[i].left = in_bit0[i];
    cmps[i].right = in_bit1[i];
  }

  // 计算所有比特中的AND运算结果。
  signal and_result = andMany(cmps);

  // 如果and结果为1，则表示in[0] 严格小于 in[1]，否则为0。
  component out_cmp = open(Cmp());
  out_cmp.left = constOfValue(1);
  out_cmp.right = and_result;
  out = out_cmp.output;
}

// 测试 LessThan 电路是否正确工作。
component testLessThan() {
  signal input[2];
  signal expectedOutput;
  signal actualOutput;

  input[0] = witness(10);
  input[1] = witness(20);
  expectedOutput = witness(1);

  actualOutput = LessThan()(input).out;

  component assert = AssertEqual(2)();
  assert.left = expectedOutput;
  assert.right = actualOutput;

  return actualOutput;
}

```

扩展1：

```
include "utils.circom";

template LessThan(k) {
  signal in[2];
  signal out;

  // 计算每个输入信号的比特数(向上取整到 k / 2)。
  constant num_bits = ceil(k / 2);
  
  // 将输入信号拆分为单独的比特(每个不超过num_bits)。
  signal in_bits[2][num_bits];
  for (i = 0; i < 2; i++) {
    for (j = 0; j < num_bits; j++) {
      in_bits[i][j] = bits(in[i], min(j+1, num_bits));
    }
  }

  // 比较并组合每个bit slice。
  component cmps[num_bits] = [0; num_bits];
  signal gt_results[num_bits];
  for (i = 0; i < num_bits; i++) {
    cmps[i] = GT();
    cmps[i].left = in_bits[0][i];
    cmps[i].right = in_bits[1][i];
    gt_results[i] = cmps[i].output;
  }

  // 使用多层 mux 组合 bit slices。
  signal final_result = muxn([1;num_bits], gt_results);
  
  // 如果final_result为1，则表示in[0] 严格小于 in[1]，否则为0。
  component out_cmp = open(Cmp());
  out_cmp.left = constOfValue(1);
  out_cmp.right = final_result;
  out = out_cmp.output;
}

// 测试 LessThan 电路是否正确工作。
component testLessThan() {
  signal input[2];
  signal expectedOutput;
  signal actualOutput;

  input[0] = witness(10);
  input[1] = witness(20);
  expectedOutput = witness(1);

  actualOutput = LessThan(8)()(input).out;

  component assert = AssertEqual(2)();
  assert.left = expectedOutput;
  assert.right = actualOutput;

  return actualOutput;
}

```

扩展2：

### 整数除法 IntegerDivide

注意：这个电路非常难！

* 参数：`nbits`。 使用 `assert` 断言这最多为 126！
* 输入信号：`dividend`, `divisor` （被除数，除数）
* 输出信号：`remainder`, `quotient` （余数，商）

要求：首先，检查 `dividend`和 `divisor`是否最多为 `nbits`位长。 接下来，计算并约束 `余数`和 `商`。

```
template IntegerDivide(nbits) {
  signal dividend, divisor;
  signal remainder, quotient;

  // 确保输入数量不超过 nbits。
  assert(dividend.numBits <= nbits, "Error: dividend must be at most nbits bits long");
  assert(divisor.numBits <= nbits, "Error: divisor must be at most nbits bits long");

  // 将被除数和除数重新定型为有 nbits 位的位数组。
  const num_chunks = ceil(nbits / 64);
  dividend_bits[num_chunks];
  divisor_bits[num_chunks];

  for (i = 0; i < num_chunks; i++) {
    dividend_bits[i] = bits(dividend, min(64, nbits - 64*i));
    divisor_bits[i] = bits(divisor, min(64, nbits - 64*i));

    component zero_extend_dividend = open(ZeroExtend(63))();
    zero_extend_dividend.input[0] = dividend_bits[i];
    dividend_bits[i] = zero_extend_dividend.output;

    component zero_extend_divisor = open(ZeroExtend(63))();
    zero_extend_divisor.input[0] = divisor_bits[i];
    divisor_bits[i] = zero_extend_divisor.output;
  }

  // 将除法器应用于两个数字中的每个比特，并使用 schoolbook 算法计算商和余数。
  // 学校算法的参考：https://en.wikipedia.org/wiki/Long_division
  const num_total_bits = num_chunks * 64;
  quotient_bits[num_total_bits];
  remainder_bits[num_total_bits];
  remainder_bits[0] = bits(dividend, min(64, nbits));
  for (i = 0; i < num_total_bits; i++) {
    cmp = GreaterEqualThan(64)()(remainder_bits[i], divisor_bits[num_chunks-1]);
    for (j = num_chunks - 2; j >= 0; j--) {
      cmp_next = GreaterEqualThan(64)()();
      component left_shift = open(LeftShift(64))();
      left_shift.input = [remainder_bits[j+1], cmp.output];
      cmp_next.in[0] = left_shift.output;
      cmp_next.in[1] = divisor_bits[j];
      cmp = And()()(cmp, cmp_next);
      component subtract = open(Subtract(64))();
      subtract.left = remainder_bits[j+1];
      subtract.right = select(cmp.output, 0)(divisor_bits[j], E()());
      remainder_bits[j+1] = subtract.output;
    }
    quotient_bits[i] = cmp.output;
    remainder_bits[0] = select(cmp.output, 0)(remainder_bits[0], E()());
    component left_shift = open(LeftShift(64))();
    left_shift.input = [remainder_bits[0], cmp.output];
    remainder_bits[0] = left_shift.output;
  }

  // 将商和余数重新合并为单个数字。
  quotient = pack(quotient_bits);
  remainder = pack(remainder_bits);

  return { remainder, quotient };
}

```

* **扩展** ：您将如何修改电路以处理负的被除数？要处理负数的被除数，我们可以采用符号位扩展。首先，我们需要将商品和余数的最高位分离出来，然后对于被除数进行符号位扩展以保持正确的值。此外，在根据符号比较商和除数时，我们还需要考虑它们的符号，因此我们需要在原来的电路的基础上添加一些逻辑门。


  * 在获取被除数和除数的符号时，我们将其与x=0作比较并使用减法处理。
  * 我们需要为底数为0的情况单独添加逻辑。
  * 最后，我们根据被除数和除数的符号确定商和余数的符号

[解决方案](https://github.com/darkforest-eth/circuits/blob/master/perlin/perlin.circom#L44)（忽略第二个参数SQRT_P，这是无关紧要的）
