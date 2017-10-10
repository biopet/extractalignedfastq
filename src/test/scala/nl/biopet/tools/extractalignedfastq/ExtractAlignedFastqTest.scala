package nl.biopet.tools.extractalignedfastq

import nl.biopet.test.BiopetTest
import org.testng.annotations.Test

object ExtractAlignedFastqTest extends BiopetTest {
  @Test
  def testNoArgs(): Unit = {
    intercept[IllegalArgumentException] {
      ToolTemplate.main(Array())
    }
  }
}
