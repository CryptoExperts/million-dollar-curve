
$(function() {
$('a[href*=#]:not([href=#])').click(function() {
  var target = $(this.hash);
  target = target.length ? target : $('[id=' + this.hash.slice(1) + ']');
  if (target.length) {
    $('html,body').animate({
      scrollTop: target.offset().top
    }, 300);
    return false;
  }
  return true;
});



$(document).ready(function() {
  var $sidebar = $('#sidebar');
  var sidebarPosition;
  var affix = function() {
    if ($(window).scrollTop() > sidebarPosition) {
      $sidebar.addClass('fixed');
    } else {
      $sidebar.removeClass('fixed');
    }
  };
  var position = function() {
    sidebarPosition = $sidebar.parent().offset().top;
    affix();
  };
  $(window).resize(position);
  $(window).scroll(affix);
  position();
  $('[data-mdc-mt]').each(function() {
    var s = $(this).attr('data-mdc-mt') + '@cryptoexperts.com';
    $(this).attr('href', 'mailto:' + s);
    if ($(this).attr('data-mdc-rt')) {
      $(this).text(s);
    }
  });
});
});

// $('#sidebar').scrollspy();
// $('#sidebar').on('activate.bs.scrollspy', function () {
//   console.log('toto');
// });

$(window).on('activate.bs.scrollspy', function () {
  var $active = $('#sidebar a.active');
  var $subnav = $active.parents('ul.subnav');
  $('#sidebar .open-active').removeClass('open-active');
  $active.parent().addClass('open-active');
  $('#sidebar a[href=' + $subnav.attr('subnavfor') + ']').parent().addClass('open-active');
});
