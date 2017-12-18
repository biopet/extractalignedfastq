package nl.biopet.tools.extractalignedfastq

import java.io.File

import nl.biopet.utils.tool.AbstractOptParser

class ArgsParser(toolCommand: ToolCommand[Args])
    extends AbstractOptParser[Args](toolCommand) {

  head(s"""
          |$cmdName - Select aligned FASTQ records
      """.stripMargin)

  opt[File]('I', "input_file") required () valueName "<bam>" action { (x, c) =>
    c.copy(inputBam = x)
  } validate { x =>
    if (x.exists) success else failure("Input BAM file not found")
  } text "Input BAM file"

  opt[String]('r', "interval") required () unbounded () valueName "<interval>" action {
    (x, c) =>
      // yes, we are appending and yes it's O(n) ~ preserving order is more important than speed here
      c.copy(intervals = c.intervals :+ x)
  } text "Interval strings (e.g. chr1:1-100)"

  opt[File]('i', "in1") required () valueName "<fastq>" action { (x, c) =>
    c.copy(inputFastq1 = x)
  } validate { x =>
    if (x.exists) success else failure("Input FASTQ file 1 not found")
  } text "Input FASTQ file 1"

  opt[File]('j', "in2") optional () valueName "<fastq>" action { (x, c) =>
    c.copy(inputFastq2 = Option(x))
  } validate { x =>
    if (x.exists) success else failure("Input FASTQ file 2 not found")
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
