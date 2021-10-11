function reddenPage() {
  var x = document.getElementsByTagName("button");
  var y = Array.prototype.filter.call(x, function(z) {
    return z.title === "Add to Card";
  });
  if(y.length > 0) {
    y[0].click();
    setTimeout(reddenPage, 2000);
  }
}

/*
function reddenPage() {
  var x = document.getElementById("myBtn");
  x.click();
  var y = document.getElementsByClassName("close");
  if(y.length > 0) {
    y[0].click();
  }
}
*/

chrome.action.onClicked.addListener((tab) => {
  chrome.scripting.executeScript({
    target: { tabId: tab.id },
    function: reddenPage
  });
});
