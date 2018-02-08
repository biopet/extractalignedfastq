/*
 * Copyright (c) 2014 Biopet
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package nl.biopet.tools.extractalignedfastq

import java.io.File

import nl.biopet.utils.tool.{AbstractOptParser, ToolCommand}

class ArgsParser(toolCommand: ToolCommand[Args])
    extends AbstractOptParser[Args](toolCommand) {

  opt[File]('I', "input_file") required () valueName "<bam>" action { (x, c) =>
    c.copy(inputBam = x)
  } text "Input BAM file"

  opt[String]('r', "interval") required () unbounded () valueName "<interval>" action {
    (x, c) =>
      // yes, we are appending and yes it's O(n) ~ preserving order is more important than speed here
      c.copy(intervals = c.intervals :+ x)
  } text "Interval strings (e.g. chr1:1-100)"

  opt[File]('i', "in1") required () valueName "<fastq>" action { (x, c) =>
    c.copy(inputFastq1 = x)
  } text "Input FASTQ file 1"

  opt[File]('j', "in2") optional () valueName "<fastq>" action { (x, c) =>
    c.copy(inputFastq2 = Option(x))
  } text "Input FASTQ file 2 (default: none)"

  opt[File]('o', "out1") required () valueName "<fastq>" action { (x, c) =>
    c.copy(outputFastq1 = x)
  } text "Output FASTQ file 1"

  opt[File]('p', "out2") optional () valueName "<fastq>" action { (x, c) =>
    c.copy(outputFastq2 = Option(x))
  } text "Output FASTQ file 2 (default: none)"

  opt[Int]('Q', "min_mapq") optional () action { (x, c) =>
    c.copy(minMapQ = x)
  } text "Minimum MAPQ of reads in target region to remove (default: 0)"

  opt[Int]('s', "read_suffix_length") optional () action { (x, c) =>
    c.copy(commonSuffixLength = x)
  } text
    """Length of suffix mark from each read pair (default: 0). This is used for distinguishing read pairs with
         different suffices. For example, if your FASTQ records end with `/1` for the first pair and `/2` for the
         second pair, the value of `read_suffix_length` should be 2."
      """.stripMargin

  note(
    """
         |This tool creates FASTQ file(s) containing reads mapped to the given alignment intervals.
       """.stripMargin)

  checkConfig { c =>
    if (c.inputFastq2.isDefined && c.outputFastq2.isEmpty)
      failure("Missing output FASTQ file 2")
    else if (c.inputFastq2.isEmpty && c.outputFastq2.isDefined)
      failure("Missing input FASTQ file 2")
    else
      success
  }
}
